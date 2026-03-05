
cat("Nested vRF Processing Log\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("Working directory:", getwd(), "\n")
cat("=", rep("=", 70), "\n\n", sep = "")

# Source utility/subfunctions (low level)
source("vRF/subfunctions_vRF.R") ## from https://github.com/mematt/ml4calibrated450k/vRF/subfunctions_vRF.R

# Source train function (mid level)
source("vRF/train_vRF.R")  ## from https://github.com/mematt/ml4calibrated450k/vRF/train_vRF.R

# Source nestedcv function (high level)
source("vRF/nestedcv_vRF.R")  ## from https://github.com/mematt/ml4calibrated450k/vRF/nestedcv_vRF.R

## Run
### Define parallel backend using `doMC` 
# Register parallel backend 
# 1. doMC  
library(doMC)
cores <- 16
# Record parallel configuration
cat(paste0("   📌 Parallel Backend: doMC\n"))
cat(paste0("   📌 Registered Cores: ", cores, "\n"))

# Test parallel backend if it is running/functioning - by Hadley Wickham - AdvR p.374. 
my_pause <- function(i){
  function(x) Sys.sleep(i)
}

# Sequential test
cat("   Running sequential test (10 x 0.25s sleep)...\n")
system.time(lapply(1:10, my_pause(0.25)))
#   user  system elapsed 
#  0.008   0.000   2.502
system.time(mclapply(1:10, my_pause(0.25), mc.cores = cores)) 
#   user  system elapsed 
#  0.004   0.272   0.269

# --------------------------
# Load data objects & record details
# --------------------------
load("data/y.RData")
load("data/nfolds.RData")
path.betas <- "/media/user/data/Sanbo_Methylation/Chordoma/GSE140686/score/CV/data/betasKk/"
betas.p <- "betaKk"
out_path <- "vRF/CVfold"
out_fname <- "CVfold"
n.cv.folds <- 3
ntrees <- 10000
p <- 10000
subset.CpGs.1k <- FALSE

# --------------------------
# Record run parameters
# --------------------------
cat("📋 Core run parameters for run_nestedcv_vRF...\n")
cat(paste0("   ➡️ path.betas.var.filtered: ", path.betas, "\n"))
cat(paste0("   ➡️ fname.betas.p.varfilt: ", betas.p, "\n"))
cat(paste0("   ➡️ n.cv.folds: ", n.cv.folds, "\n"))
cat(paste0("   ➡️ out.path: ", out_path, "\n"))
cat(paste0("   ➡️ out.fname: ", out_fname, "\n"))
cat(paste0("   ➡️ subset.CpGs.1k: ", subset.CpGs.1k, "\n"))
cat(paste0("   ➡️ cores: ", cores, "\n"))
cat(paste0("   ➡️ ntrees: ", ntrees, "\n"))
cat(paste0("   ➡️ p: ", p, "\n"))
cat("\n")

# --------------------------
# Execute core function & record timing
# --------------------------
cat("Starting run_nestedcv_vRF execution...\n")
start_time <- Sys.time()
cat("   Execution start time:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")

run_nestedcv_vRF(
    path.betas.var.filtered = path.betas,
    fname.betas.p.varfilt = betas.p,
    n.cv.folds = n.cv.folds,
    ntrees = ntrees, p = p,
    K.start = 1, k.start = 0,
    out.path = out_path, out.fname = out_fname,
    subset.CpGs.1k = subset.CpGs.1k,
    cores = cores)
end_time <- Sys.time()
total_duration <- difftime(end_time, start_time, units = "mins")  
cat("Finished run_nestedcv_vRF execution.\n")
cat("Execution end time:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")


message("Log file saved to: ", normalizePath(log_file))
