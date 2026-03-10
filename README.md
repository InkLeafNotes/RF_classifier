# RF_classifier
To analyze the DNA methylation profile in the context of sarcoma classifier,  random forest (RF)-based sarcoma classification model based on the data of the sarcoma classifier reference cases was reconstructed, this enabled the generation of calibrated scores for individual samples.

This repository contains R-scripts used to perform DNA-methylation data analysis on pediatric chordoma. The raw data of sarcoma classifier reference set is publicly available at  GEO under Accession number [GSE140686](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140686).

## Preprocessing reference datset
Reference data preprocessing (normalization, probe filtering, batch adjustment, annotation merging)
```R
preprocess.R
```

## Training Classifier 
Train the RF classification model with parallel computing and feature selection
```R
RF_model/run_train_RF_model.R
```

## Cross Validation
### Data Preprocessing
Process reference dataset into nested CV-compatible data structure (fold generation, data splitting/matching)
```R
data/run_makebetaKk.R
```

### Internal validation
Nested cross-validation for RF model performance evaluation
```R
vRF/run_nested_vRF.R
```

## MR model Calibration
Train and save the final universal multinomial ridge regression (MR) model for score calibration
```R
calibrate/calibration_MR_model.R
```
Run MR model calibration workflow, the results can be used for the performance_evaluator in ml4calibrated450k
```R
calibrate/run_calibrate_vRF_MR.R
```

## Predict New Samples
Preprocess Test Data
850K methylation array data preprocessing (normalization, batch effect correction, beta value calculation)
```R
predict/preprocess.R
```
Predict with trained RF model and calibrate scores using the MR model, output integration results
```R
predict/predict.R
```
## Key Dependencies
minfi, randomForest, glmnet, doMC, data.table, GEOquery, limma

## Third-party Code Reference
Scripts adapted from [ml4calibrated450k](https://github.com/mematt/ml4calibrated450k) and [mnp_training](https://github.com/mwsill/mnp_training) (see individual script comments for details)
