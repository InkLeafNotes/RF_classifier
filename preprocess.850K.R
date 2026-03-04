setwd("/media/user/data/Sanbo_Methylation/Chordoma/GSE140686/score/predict")
library(minfi)

source("MNPprocessIDAT_functions.R")

##loading data
baseDir <- "."
targets <- read.metharray.sheet(baseDir)
idatdir <- "/media/user/data/Sanbo_Methylation/rawdata/"
targets$Basename <- paste(idatdir,"idat/",targets$Slide,"_",targets$Array,sep = "")

rgSet <- read.metharray.exp(targets=targets,force=T)
sampleNames(rgSet) <- targets$Sample_Name

## Illumina normalization
message("running normalization ...",Sys.time())
Mset <- MNPpreprocessIllumina(rgSet)

## get betas directly
methy <- getMeth(Mset)
unmethy <- getUnmeth(Mset)
betas <- methy / (methy +unmethy +100) # same as getBeta(Mset, offset = 100)
betas.nmp  <- as.data.frame(t(betas))
saveRDS(betas.nmp, "chordoma.nmp.rds")

## batch correction and then get betas
library(limma)
message("running batch correction ...",Sys.time())
load("/media/user/data/GSE/GSE140686/preprocess/results/Mset_filtered.RData")
# get FFPE/Frozen type
ffpe <- anno$`material preparation:ch1`
ffpe <- c(ffpe,rep("FFPE",length(targets$Sample_Name)))
batch <- ifelse(ffpe == "FFPE", 2, 1)

# combine Mset and Mset_filtered
# 1. get probe ID of two Mset
probes_Mset <- rownames(Mset)
probes_Mset_filtered <- rownames(Mset_filtered)

# 2. Filter for common probes (intersection)
common_probes <- intersect(probes_Mset, probes_Mset_filtered)
cat("两个 Mset 的共有探针数量：", length(common_probes), "\n")

# 3. Filter both Msets to retain only common probes 
Mset <- Mset[common_probes, ]
Mset_filtered <- Mset_filtered[common_probes, ]

# 4. combined
mset_combined <- combine(Mset_filtered, Mset)

# Verify the combination result
cat("Number of samples in Mset_filtered before combination: ", ncol(Mset_filtered), "\n")
cat("Number of probes in Mset_filtered after combination: ", nrow(Mset_filtered), "\n")
cat("Number of samples in Mset before combination: ", ncol(Mset), "\n")
cat("Number of probes in Mset after combination: ", nrow(Mset), "\n")
cat("Number of samples in mset_combined after combination: ", ncol(mset_combined), "\n")
cat("Number of probes in mset_combined after combination: ", nrow(mset_combined), "\n")

methy_c <- getMeth(mset_combined)
unmethy_c <- getUnmeth(mset_combined)
# remove batch effects by linear model
methy.ba <- 2^removeBatchEffect(log2(methy_c +1), batch)
unmethy.ba <- 2^removeBatchEffect(log2(unmethy_c +1), batch)

# recalculate betas, illumina like
betas_c <- methy.ba / (methy.ba +unmethy.ba +100)
chord_only <- setdiff(colnames(betas_c), colnames(Mset_filtered))
betas_ch <- as.data.frame(t(betas_c[, chord_only]))
saveRDS(betas_ch, "chordoma.ba.rds")
message("preprocessing finished ...",Sys.time())


## 用reference batch coeficients进行批次处理 

load("/media/user/data/GSE/GSE140686/preprocess/results/ba.coef.RData")

# 1. get probe ID 
probes_Mset <- rownames(Mset)
probes_reference <- names(methy.coef$FFPE)

# 2. Filter for common probes (intersection)
common_probes <- intersect(probes_Mset, probes_reference )

# 3. Filter both Msets to retain only common probes 
Mset <- Mset[common_probes, ]
cat(" Mset和reference共有探针数量：", nrow(Mset), "\n")

batch <- rep("FFPE",length(Mset$Sample_Name))

methy <- getMeth(Mset)
unmethy <- getUnmeth(Mset)

# perform batch adjustment
methy.b <- log2(methy +1) + matrix(unlist(methy.coef[match(batch,names(methy.coef))]),ncol=length(batch))
unmethy.b <- log2(unmethy +1) + matrix(unlist(unmethy.coef[match(batch,names(unmethy.coef))]),ncol=length(batch))
methy.b[methy.b < 0] <- 0
unmethy.b[unmethy.b < 0] <- 0
methy.ba <- 2^methy.b
unmethy.ba <- 2^unmethy.b
# illumina-like beta values
betas.test <- methy.ba / (methy.ba +unmethy.ba +100)
betas.test <- as.data.frame(t(betas.test))
saveRDS(betas.test, "chordoma.ba.coef.rds")
