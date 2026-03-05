
# Source subfunctions for multinomial ridge regression (MR)
source("calibrate/calibration_MR.R") ## from https://github.com/mematt/ml4calibrated450k/blob/master/calibrate/calibration_MR.R
load("data/y.RData")
load("data/nfolds.RData")

# Define parameters
out.path <- "/calibrate/vRF-calibrated-MR/"
out.fname <- "probsCVfold"
load.path.w.name <- "./vRF/CVfold/CVfold."
verbose.messages <- TRUE
which.optimized.metric.or.algorithm <- "vanilla"
save.metric.name.into.output.file <- TRUE
parallel.cv.glmnet <- TRUE
setseed <- 1234

# Log parameters
cat("Parameters used for calibrate_MR:\n")
cat("===============================\n")
cat("out.path:", out.path, "\n")
cat("out.fname:", out.fname, "\n")
cat("load.path.w.name:", load.path.w.name, "\n")
cat("verbose.messages:", verbose.messages, "\n")
cat("which.optimized.metric.or.algorithm:", which.optimized.metric.or.algorithm, "\n")
cat("save.metric.name.into.output.file:", save.metric.name.into.output.file, "\n")
cat("parallel.cv.glmnet:", parallel.cv.glmnet, "\n")
cat("setseed:", setseed, "\n")
cat("===============================\n\n")

start_time <- Sys.time()
cat("Execution start time:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
calibrate_MR(out.path = out.path,
             out.fname = out.fname,
             load.path.w.name = load.path.w.name,
             verbose.messages = verbose.messages,
             which.optimized.metric.or.algorithm = which.optimized.metric.or.algorithm, 
             save.metric.name.into.output.file = save.metric.name.into.output.file,
             parallel.cv.glmnet = parallel.cv.glmnet,
             setseed = setseed)
end_time <- Sys.time()
total_duration <- difftime(end_time, start_time, units = "mins")  
cat("Finished run_nestedcv_vRF execution.\n")
cat("Execution end time:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")


