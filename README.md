# RF_classifier
To analyze the DNA methylation profile in the context of sarcoma classifier,  random forest (RF)-based sarcoma classification model based on the data of the sarcoma classifier reference cases was reconstructed, this enabled the generation of calibrated scores for individual samples.

Collection of R-scripts used to perform DNA-methylation data analysis on pediatric chordoma. The raw data of sarcoma classifier reference set is publicly available at  GEO under Accession number [GSE140686](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140686).

## Preprocessing reference datset
```R
preprocess.R
```

## Training Classifier 
```R
RF_model/run_train_RF_model.R
```

## Cross Validation
### Data Preprocessing
```R
data/run_makebetaKk.R
```

### Internal validation using RF classifier
```R
vRF/run_nested_vRF.R
```

## MR model Calibration
```R
calibrate/calibration_MR_model.R
```

## Predict
Preprocess the 850K methylation array data
```R
predict/preprocess.R
```
predict new sample calibrated score with RF-model and MR-model
```R
predict.R
```
