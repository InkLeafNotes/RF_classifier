## Run
### Define parallel backend using `doMC` 
# Register parallel backend --------------------------------------------------------------------------------------------------------------------------------------
# 1. doMC  ----------------------------------------------------------------------------------------------------------------------
library(doMC)
cores <- 10
# Record parallel configuration
cat(paste0("   📌 Parallel Backend: doMC\n"))
cat(paste0("   📌 Registered Cores: ", cores, "\n"))
# Test parallel backend if it is running/functioning - by Hadley Wickham - AdvR p.374. -------------------------------------------
my_pause <- function(i){
  function(x) Sys.sleep(i)
}

# Sequential test
cat("   Running sequential test (10 x 0.25s sleep)...\n")
system.time(lapply(1:10, my_pause(0.25)))
#   user  system elapsed 
#  0.008   0.000   2.502
system.time(mclapply(1:10, my_pause(0.25), mc.cores = cores)) # 10
#   user  system elapsed 
#  0.004   0.272   0.269

# --------------------------
# Load data objects & record details
# --------------------------

library(randomForest)
source("train_vRF.R") ## from Github(https://github.com/mematt/ml4calibrated450k)

cat("📊 Starting input data preparation...\n")
load("/media/user/data/GSE/GSE140686/preprocess/results/betas.ba.RData") ## including betas(reference beta value matrix) and anno(sample anno)

cat(paste0("✅ Reference beta data loaded (time taken: ", round(difftime(end_read, start_read, units = "secs"), 2), "s)\n"))
cat(paste0("   - Reference beta dimensions: ", nrow(betas), " samples × ", ncol(betas), " probes\n"))


cat("\n🧹 Filtering probes with NA values...\n")
na_per_probe <- colSums(is.na(betas))
valid_probes <- colnames(betas)[colSums(is.na(betas)) == 0]
betas <- betas[, valid_probes, drop = FALSE]
cat(paste0("⚙️ Total number of probes used for feature selection: ", ncol(betas)))

anno$Group <- as.factor(anno$Group)

cat(paste0("   - Sample info dimensions: ", nrow(anno), " samples × ", ncol(anno), " columns\n"))
cat(paste0("   - Number of unique groups: ", length(unique(anno$Group)), "\n"))

cat("\n🔍 Validating sample matching between beta data and sample info...\n")
sample_match <- identical(anno$PID, rownames(betas))
cat(paste0("   - All sample IDs in anno exist in beta data: ", ifelse(sample_match, "✅ YES", "❌ NO"), "\n"))

if (!sample_match) {
  mismatched_samples <- setdiff(anno$PID, row.names(betas))
  cat(paste0("   - Mismatched sample IDs: ", paste(mismatched_samples, collapse = ", "), "\n"))
}

ntrees <- 10000
p <- 10000
seed <- 1234

# --------------------------
# Record run parameters 
# --------------------------
cat("📋 Core run parameters for run_train_RF...\n")
cat(paste0("   ➡️ cores: ", cores, "\n"))
cat(paste0("   ➡️ ntrees: ", ntrees, "\n"))
cat(paste0("   ➡️ p: ", p, "\n"))
cat("\n")

# train RF
# --------------------------
# Execute core function & record timing
# --------------------------
cat("Starting run_train_RF execution...\n")
start_time <- Sys.time()
cat("   Execution start time:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")

rfcv <- trainRF(y = anno$Group,
                      betas = betas, 
                      ntrees = ntrees,
                      p = p,
                      seed = seed,
                      cores = cores)

# Save RF 
saveRDS(rfcv, file = paste0("RF_model", format(Sys.time(), "%Y%m%d_%H%M"), ".rds"))

cat("Finished run_train_RF execution.\n")

