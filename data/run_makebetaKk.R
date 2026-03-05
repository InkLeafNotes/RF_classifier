#-------------------------------------------------
# Process the reference dataset into a nested data structure suitable for cross-validation (CV).
#
# 2026-01-12
#-------------------------------------------------

cat("🔧 Loading required scripts...\n")
source("./data/makefolds.R") ## from https://github.com/mematt/ml4calibrated450k/data/makefolds.R
source("./data/load_match_save_betaKk.R")
cat("✅ Scripts loaded successfully\n\n")

# --------------------------
# load reference data
# --------------------------
cat("📊 Starting input data preparation...\n")
load("/media/user/data/GSE/GSE140686/preprocess/results/betas.ba.RData") ## including betas(reference beta value matrix) and anno(sample anno)

cat(paste0("   - Reference beta dimensions: ", nrow(betas), " samples × ", ncol(betas), " probes\n"))

anno$Group <- as.factor(anno$`methylation class:ch1`)

cat(paste0("   - Sample info dimensions: ", nrow(anno), " samples × ", ncol(anno), " columns\n"))
cat(paste0("   - Number of unique groups: ", length(unique(anno$Group)), "\n"))

cat("\n🔍 Validating sample matching between beta data and sample info...\n")
sample_match <- identical(anno$geo_accession, rownames(betas))
cat(paste0("   - All sample IDs in anno exist in beta data: ", ifelse(sample_match, "✅ YES", "❌ NO"), "\n"))

if (!sample_match) {
  mismatched_samples <- setdiff(anno$geo_accession, row.names(betas))
  cat(paste0("   - Mismatched sample IDs: ", paste(mismatched_samples, collapse = ", "), "\n"))
}


anno <- as.data.frame(anno, row.names = as.character(anno$geo_accession))
save(anno, file = "data/anno.RData")
y <- anno$Group
save(y, file = "data/y.RData")
cat("✅ Sample annotation saved to data/anno.RData and data/y.RData\n")


cat("\n🧹 Filtering probes with NA values...\n")
na_per_probe <- colSums(is.na(betas))
valid_probes <- colnames(betas)[colSums(is.na(betas)) == 0]
cat(paste0("   - Original number of probes: ", ncol(betas), "\n"))
cat(paste0("   - Valid probes (no NA): ", length(valid_probes), "\n"))

betas <- betas[, valid_probes, drop = FALSE]
cat(paste0("⚙️ Total number of probes used for feature selection: ", ncol(betas)))

# Inner/Nested 3-fold CV loops
cat("\n🔄 Generating nested 3-fold CV folds...\n")
seed <- 1234
cat("setseed:",seed,"\n")
set.seed(seed)
nfolds <- makenestedfolds(anno$Group, cv.fold = 3)
save(nfolds, file = "data/nfolds.RData")
cat("✅ Nested CV folds generated and saved to data/nfolds.RData\n")

betas.mat <- t(betas)
cat("Matrix dimensions:", dim(betas.mat)[1], "CpGs ×", dim(betas.mat)[2], "samples\n")

load_match_save_betasKk(
    betas.mat = betas.mat,
    anno = anno,
    nfolds.. = nfolds,
    n.cv.folds = 3,
    subset.10K = FALSE,
    out.path = "data/betasKk",
    out.fname = "betaKk")

# Function call example
# Sys.time()
# load_filter_match_save_betasKk(
#   K.start = 1,
#   k.start = 0,
#   n.cv.folds = 5,
#   nfolds.. = your_cv_folds_list,  # Your cross-validation folds object
#   betas.mat = your_beta_matrix,   # Your beta matrix (rows=CpGs, columns=samples)
#   anno = your_annotation_data,    # Annotation data with sample IDs in $X column
#   out.path = "betas.varfilt.10k", # Output directory
#   out.fname = "betas.K.k",         # Output filename prefix
#   subset.10K = FALSE             # Subset to top 10K probes to test
# )
# Sys.time()

