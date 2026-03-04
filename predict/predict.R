# --------------------------
# 设置日志
# --------------------------
log_dir <- "logs"
dir.create(log_dir, showWarnings = FALSE)
log_file <- file.path(log_dir, paste0("train_RF_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))

# 开始记录
message("Starting nested CV processing. Log will be saved to: ", log_file)
con <- file(log_file, open = "wt")
sink(con, type = "output")
sink(con, type = "message")

library(randomForest)
## load RF model
rfcv <- readRDS("../RF_model/RF_model.rds")

## load test data
# 行是样本，列是探针
beta.test <- readRDS("chordoma.nmp.rds")

beta.test <- readRDS("chordoma.ba.rds")

beta.test <- readRDS("chordoma.ba.coef.rds")


# 1. 提取模型的特征名（用于匹配）
model_feature_names <- rownames(rfcv[[1]]$importance)

# 2. 校验特征名匹配情况（关键！）
col_matches <- match(model_feature_names, colnames(beta.test))

# 3. 过滤掉匹配失败的列（避免NA索引）
valid_col_idx <- col_matches[!is.na(col_matches)]
valid_feature_names <- model_feature_names[!is.na(col_matches)]

if (length(valid_col_idx) == 0) {
  stop("❌ No features matched between model and test data!")
} else {
  cat(sprintf("\n✅ Filtered to %d valid matched features (out of %d)\n", 
              length(valid_col_idx), length(model_feature_names)))
}

# 4. 提取测试数据并强制保持数据框格式（关键！）
test_data <- beta.test[, valid_col_idx, drop = FALSE]

# 5. 二次检查测试数据的缺失值（包含隐性问题）
if (sum(is.na(test_data)) > 0) {
  cat("\n⚠️ Found NA values in test data, imputing with median...\n")
  # 用中位数填充（避免删除样本）
  test_data <- apply(test_data, 2, function(x) {
    x[is.na(x)] <- median(x, na.rm = TRUE)
    return(x)
  })
  test_data <- as.data.frame(test_data)
  colnames(test_data) <- valid_feature_names
}

# 6. 执行预测（最终修复版）
cat("\n🚀 Starting prediction with clean test data...\n")
scores <- predict(rfcv[[1]],
                  test_data,
                  type = "prob")

cat("\n✅ Prediction completed successfully! First 5 rows of scores:\n")
#print(head(scores))

# 创建结果表格
# 1. 提取每个样本的最终分类标签（基于最大概率）
predicted_class <- colnames(scores)[max.col(scores)]

# 2. 提取最终类别对应的原始分数
raw_score_pred_class <- sapply(1:nrow(scores), function(i) {
  scores[i, predicted_class[i]]  # 第i行，取“最终预测类别”列的数值
})

# 3. 将分类结果和概率合并成表格（方便查看）
score_with_class <- data.frame(
  样本ID = rownames(scores),  # 保留样本ID（如果有）
  预测分类 = predicted_class,
  概率 = raw_score_pred_class,
  scores  # 保留原始概率值
)

# 4. 查看合并后的结果
cat("\n📋 带分类标签的完整结果（前5行）：\n")
#print(head(score_with_class))

library(glmnet)
# 关键：将原始分数转换为矩阵（glmnet要求输入为矩阵）
rf_raw_probs_mat <- as.matrix(scores)
load("../Archives/calibrate/final_mr_calibrator.RData")
load("../CV/calibrate/final_mr_calibrator.RData")
# 用MR模型校正分数（核心步骤）
calibrated_probs <- predict(
  object = cv_final_mr$glmnet.fit,  # MR模型的核心对象（不是cv_final_mr本身！）
  newx = rf_raw_probs_mat,          # 新样本的原始RF分数（矩阵）
  type = "response",                # 输出校正后的概率（不是系数）
  s = cv_final_mr$lambda.1se        # 最优lambda（用1se更稳健，也可换lambda.min）
)[,,1]  # [,,1] 消除glmnet输出的三维数组维度（仅保留样本×类别矩阵）

# 核心修复：为校正后的概率矩阵添加列名
colnames(calibrated_probs) <- colnames(rf_raw_probs_mat) # 继承原始RF的类别列名

# ========== 第四步：结果整合（同之前） ==========
# 生成最终预测类别（此时列名存在，可正确匹配）
final_pred_classes <- colnames(calibrated_probs)[apply(calibrated_probs, 1, which.max)]

# 1. 提取最终类别对应的原始分数
raw_score <- sapply(1:nrow(rf_raw_probs_mat), function(i) {
  rf_raw_probs_mat[i, final_pred_classes[i]]  # 第i行，取“最终预测类别”列的数值
})

# 2. 提取最终类别对应的校正分数
calib_score <- sapply(1:nrow(calibrated_probs), function(i) {
  calibrated_probs[i, final_pred_classes[i]]  # 同上
})
# 整合结果
result_df <- data.frame(
  Sample_ID = rownames(scores),
  Final_Pred_Class = final_pred_classes,
  Raw_Score = raw_score,
  Calibrated_Score = calib_score,
  rf_raw_probs_mat,                       # 原始分数（带列名）
  calibrated_probs                    # 校正分数（已补列名）
)

write.table(result_df, file = "result.ba.co.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# --------------------------
# 关闭日志
# --------------------------
sink(type = "output")
sink(type = "message")
close(con)
message("Nested CV processing finished. Log saved to: ", log_file)
