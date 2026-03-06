#-----------------------------------------------------------------------------------
# preprocessing 

library(minfi)
library(GEOquery)
library(limma)
library(openxlsx)
library(tidyr)
library(dplyr)

# Set up log directory and file
log_dir="./logs"
log_file <-file.path(log_dir, paste0("preprocess_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
# Create logs directory if it doesn't exist
if (!dir.exists(log_dir)) {
  dir.create(log_dir, recursive = TRUE)
}

# Open connection and start logging
con <- file(log_file, open = "wt")
sink(con, type = "output")
sink(con, type = "message")

# Set up results directory
if (!dir.exists("results")) dir.create("results")

cat("Starting preprocessing script\n")
cat("Working directory:", getwd(), "\n")

# Source functions
source(file.path("R","MNPprocessIDAT_functions.R")) ## from https://github.com/mwsill/mnp_training/tree/master/R/MNPprocessIDAT_functions.R
cat("MNPprocessIDAT_functions.R sourced successfully\n")

# Load annotation data
table <- read.xlsx("../Table S1.xlsx") ## from https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-20603-4/MediaObjects/41467_2020_20603_MOESM4_ESM.xlsx
cat("Table S1.xlsx loaded successfully\n")

# Get sample annotation from GEO
cat("Fetching GEO data...\n")
gse <- getGEO("GSE140686", GSEMatrix=TRUE, getGPL=FALSE)
cat("GEO data fetched successfully\n")

# Process annotations
anno1 <- pData(gse[[1]])
anno2 <- pData(gse[[2]])

anno1$platform <- names(gse)[1]
anno2$platform <- names(gse)[2]

anno <- rbind(anno1, anno2)
cat("Annotations processed successfully\n")

# Load additional annotation files
ann <- read.table("../GSE140686.anno.txt", sep = "\t", header = T)
cat("GSE140686.anno.txt loaded successfully\n")

# Merge annotations
merged_ann <- merge(ann, anno, by.x = "PID", by.y = "description.1")
merged_anno <- merge(merged_ann, table, by.x = "PID", by.y = "ID")
merged_anno$id <- as.numeric(gsub("REFERENCE_SAMPLE ", "", merged_anno$PID))
merged_anno <- merged_anno[order(merged_anno$id), ]

write.table(merged_anno, "GSE140686.merged.reference.anno.txt", sep = "\t", row.names = F, quote = FALSE)
cat("Annotations merged and saved successfully\n")

# Process selected annotation columns
target_cols <- c("PID", "Methylation.Class.Name.x", "Group", "geo_accession",
                 "platform_id", "material preparation:ch1", "IDAT", "Colour")

selected_anno <- merged_anno[, target_cols]

selected_anno <- selected_anno %>%
    separate(
        col = IDAT,
        into = c("Sentrix", "Sentrix_Position"),
        sep = "_",
        remove = F
    )

selected_anno$Sentrix_ID <- paste(selected_anno$geo_accession, selected_anno$Sentrix, sep = "_")
selected_anno <- selected_anno %>% select(PID, Sentrix_ID, Sentrix_Position, everything())

# Write platform-specific sample lists
anno_p1 <- selected_anno[selected_anno$platform_id == "GPL13534", ]
write.csv(anno_p1, "sample.list.GPL13534.csv", row.names = F, quote = FALSE)

anno_p2 <- selected_anno[selected_anno$platform_id == "GPL21145", ]
write.csv(anno_p2, "sample.list.GPL21145.csv", row.names = F, quote = FALSE)

cat("Platform-specific sample lists created successfully\n")

# IDAT directory
idat_dir <- "../RAW/"
base_dir <- "."

# Process 450K data
cat("Loading 450K data...\n")
targets.450K <- read.metharray.sheet(base_dir, pattern = "GPL13534")
targets.450K$Basename <- paste(idat_dir, targets.450K$Slide, "_", targets.450K$Array, sep = "")
RGset.450K <- read.metharray.exp(targets = targets.450K, force = T)
sampleNames(RGset.450K) <- targets.450K$PID
cat("450K data loaded successfully\n")

# Process 850K data
cat("Loading 850K data...\n")
targets.850K <- read.metharray.sheet(base_dir, pattern = "GPL21145")
targets.850K$Basename <- paste(idat_dir, targets.850K$Slide, "_", targets.850K$Array, sep = "")
RGset.850K <- read.metharray.exp(targets = targets.850K, force = T)
sampleNames(RGset.850K) <- targets.850K$PID
cat("850K data loaded successfully\n")

# Combine arrays
RGset <- combineArrays(RGset.850K, RGset.450K)
cat("Arrays combined successfully\n")

# sort by id
RGset$id  <- as.numeric(gsub("REFERENCE_SAMPLE ", "", RGset$PID))
RGset <- RGset[,order(RGset$id)]

# Illumina normalization
cat("Starting normalization...\n")
Mset <- MNPpreprocessIllumina(RGset)
cat("Normalization completed successfully\n")

# Probe filtering
# all files from https://github.com/mwsill/mnp_training/tree/master/filter
cat("Starting probe filtering...\n")
amb.filter <- read.table(file.path("filter","amb_3965probes.vh20151030.txt"),header=F)
epic.filter <- read.table(file.path("filter","epicV1B2_32260probes.vh20160325.txt"),header=F)
snp.filter <- read.table(file.path("filter","snp_7998probes.vh20151030.txt"),header=F)
xy.filter <- read.table(file.path("filter","xy_11551probes.vh20151030.txt"),header=F)

rs.filter <- grep("rs", rownames(Mset))
ch.filter <- grep("ch", rownames(Mset))

# filter CpG probes
remove <- unique(c(match(amb.filter[,1], rownames(Mset)),
                   match(epic.filter[,1], rownames(Mset)),
                   match(snp.filter[,1], rownames(Mset)),
                   match(xy.filter[,1], rownames(Mset)),
                   rs.filter,
                   ch.filter))
remove <- unique(na.omit(remove))

has_na <- any(is.na(remove))
if (has_na) {
    cat("Warning: remove vector contains NA values!\n")
} else {
    cat("remove vector processed successfully - no NA values\n")
}

Mset_filtered <- Mset[-remove,]
cat("Probe filtering completed successfully\n")

# check sample matching
anno <- merged_anno
cat("\n Validating sample matching between Mset data and sample info...\n")
sample_match <- identical(anno$PID,colnames(Mset_filtered))
cat(paste0("   - Is the order of the annotations consistent with the sample data in the Mset data? : ", ifelse(sample_match, "YES", "NO"), "\n"))

# Save filtered data
save(Mset_filtered, anno, file=file.path("results","Mset_filtered.RData"))  
cat("Filtered data saved successfully\n")

# Clean up memory
rm(Mset)
gc()

# Batch adjustment
cat("Starting batch adjustment...\n")
methy <- getMeth(Mset_filtered)
unmethy <- getUnmeth(Mset_filtered)

# get FFPE/Frozen type
ffpe <- Mset_filtered$material.preparation.ch1
batch <- ifelse(ffpe == "FFPE", 2, 1)

# remove batch effects by linear model
methy.ba <- 2^removeBatchEffect(log2(methy +1), batch)
unmethy.ba <- 2^removeBatchEffect(log2(unmethy +1), batch)

# extract effects to adjust diagnostic samples
s.frozen <- min(which(batch == 1))
s.ffpe <- min(which(batch == 2))
methy.coef <- unmethy.coef <- list()
methy.coef[["Frozen"]] <- log2(methy.ba[, s.frozen]) - log2(methy[, s.frozen] +1)
methy.coef[["FFPE"]] <- log2(methy.ba[, s.ffpe]) - log2(methy[, s.ffpe] +1)
unmethy.coef[["Frozen"]] <- log2(unmethy.ba[, s.frozen]) - log2(unmethy[, s.frozen] +1)
unmethy.coef[["FFPE"]] <- log2(unmethy.ba[, s.ffpe]) - log2(unmethy[, s.ffpe] +1)

# save batch effects 
save(methy.coef,unmethy.coef,file=file.path("results","ba.coef.RData")) 

# recalculate betas, illumina like
betas <- methy.ba / (methy.ba +unmethy.ba +100)
betas <- as.data.frame(t(betas))

cat("betas dims:", dim(betas),"\n")


# check sample matching
cat("\n Validating sample matching between data and sample info...\n")
sample_match <- identical(anno$PID, rownames(betas))
cat(paste0("   - Is the order of the annotations consistent with the sample data in the betas?: ", ifelse(sample_match, "YES", "NO"), "\n"))

# save beta data
save(betas,anno,file=file.path("results","betas.ba.RData")) 



message("Log file saved to: ", normalizePath(log_file))
