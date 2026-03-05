#---------------------------------------------
# Train and save the final calibrationMR model.
# 2026-01-15
#---------------------------------------------
library(glmnet)

load("data/y.RData")
load("data/nfolds.RData")
load.path.w.name = "./vRF/CVfold/CVfold."

all_rf_scores <- list()
all_labels <- c()
  
for(i in 1:length(nfolds)){
  env_outer <- environment()
  load(paste0(load.path.w.name, i, ".", 0, ".RData"), envir = env_outer)
  scores_outer <- get("scores", envir = env_outer)
  idx_outer <- env_outer$fold$test
  all_rf_scores[[paste0("outer_",i)]] <- scores_outer
  all_labels <- c(all_labels, y[idx_outer])
}
  
rf_scores_combined <- do.call(rbind, all_rf_scores)
stopifnot(nrow(rf_scores_combined) == length(all_labels))
  
#  train final MR model
setseed = 1234
set.seed(setseed, kind = "L'Ecuyer-CMRG")
suppressWarnings({
  cv_final_mr <- cv.glmnet(
    x = rf_scores_combined,
    y = all_labels,
    family = "multinomial",
    alpha = 0,  
    type.measure = "mse",
    nlambda = 100,
    lambda.min.ratio = 1e-6,
    parallel = TRUE
  )
})

save(cv_final_mr, file = "./calibrate/final_mr_calibrator.RData")
message("最终MR模型已保存为 ./calibrate/final_mr_calibrator.RData")
