library(randomForest)
## load RF model
rfcv <- readRDS("../RF_model/RF_model.rds")

## load test data
beta.test <- readRDS("chordoma.ba.coef.rds")

model_feature_names <- rownames(rfcv[[1]]$importance)
col_matches <- match(model_feature_names, colnames(beta.test))
valid_col_idx <- col_matches[!is.na(col_matches)]
valid_feature_names <- model_feature_names[!is.na(col_matches)]
test_data <- beta.test[, valid_col_idx, drop = FALSE]

# check NA
if (sum(is.na(test_data)) > 0) {
  cat("\n⚠️ Found NA values in test data, imputing with median...\n")
  test_data <- apply(test_data, 2, function(x) {
    x[is.na(x)] <- median(x, na.rm = TRUE)
    return(x)
  })
  test_data <- as.data.frame(test_data)
  colnames(test_data) <- valid_feature_names
}

## predict
cat("\n🚀 Starting prediction with clean test data...\n")
scores <- predict(rfcv[[1]],
                  test_data,
                  type = "prob")

## calibrate scores
library(glmnet)
rf_raw_probs_mat <- as.matrix(scores)
# load calibrator Model
load("../CV/calibrate/final_mr_calibrator.RData")

calibrated_probs <- predict(
  object = cv_final_mr$glmnet.fit,  
  newx = rf_raw_probs_mat,          
  type = "response",                
  s = cv_final_mr$lambda.1se       
)[,,1]  # [,,1] 

colnames(calibrated_probs) <- colnames(rf_raw_probs_mat)

## Integration Results Table
final_pred_classes <- colnames(calibrated_probs)[apply(calibrated_probs, 1, which.max)]

raw_score <- sapply(1:nrow(rf_raw_probs_mat), function(i) {
  rf_raw_probs_mat[i, final_pred_classes[i]]  
})
calib_score <- sapply(1:nrow(calibrated_probs), function(i) {
  calibrated_probs[i, final_pred_classes[i]]  
})

result_df <- data.frame(
  Sample_ID = rownames(scores),
  Final_Pred_Class = final_pred_classes,
  Raw_Score = raw_score,
  Calibrated_Score = calib_score,
  rf_raw_probs_mat,                       
  calibrated_probs                    
)

write.table(result_df, file = "result.ba.co.txt", sep = "\t", quote = FALSE, row.names = FALSE)

