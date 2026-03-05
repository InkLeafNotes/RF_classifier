#--------------------------------------------------------------------
# Loading and variance filtering of reference betas
# modified from https://github.com/mematt/ml4calibrated450k/subfunction_load_subset_filter_match_betasKk.R
# 2026-01-12-19 UTC
#--------------------------------------------------------------------
# Function call example
# Sys.time()
# load_match_save_betasKk(
#   K.start = 1,
#   k.start = 0,
#   n.cv.folds = 5,
#   nfolds.. = your_cv_folds_list,  # Your cross-validation folds object
#   betas.mat = your_beta_matrix,   # Your beta matrix (rows=CpGs, columns=samples)
#   anno = your_annotation_data,    # Annotation data with sample IDs in $X column
#   out.path = "betas.varfilt.10k", # Output directory
#   out.fname = "betas.K.k"         # Output filename prefix
# )
# Sys.time()
#--------------------------------------------------------------------

# Check and load required packages with warning only (no installation)
if(!requireNamespace("data.table", quietly = TRUE)) {
  warning("Package 'data.table' is required but not installed.")
} else {
  library(data.table)
}


load_match_save_betasKk <- function(K.start = 1, k.start = 0, n.cv.folds = 5,
                                                              nfolds.. = NULL,
                                                              betas.mat= NULL,
                                                              anno = NULL,
                                                              subset.10K = FALSE,
                                                              out.path = "betasKk", out.fname = "betas.K.k"){

  if(!all(anno$X %in% colnames(betas.mat))) {
    warning("Some sample IDs do not exist in the matrix, only intersecting samples will be retained")
    keep_samples <- intersect(anno$X, colnames(betas.mat))
    betas.mat <- betas.mat[, keep_samples]
    anno <- anno[anno$X %in% keep_samples, ]
  }
  message("Matrix loading completed, final dimensions: rows (CpG)=", nrow(betas.mat), ", columns (samples)=", ncol(betas.mat))

 message("\nNested CV process started @ ", Sys.time())
  for(K in K.start:n.cv.folds){  
    for(k in k.start:n.cv.folds){ 
      if(k > 0){ 
        message("\n Subsetting & filtering inner/nested fold ", K,".", k,"  ... ",Sys.time())
        fold <- nfolds..[[K]][[2]][[k]]  
      } else{                                                                          
        message("\n \nSubsetting & filtering outer fold ", K,".0  ... ",Sys.time())
        fold <- nfolds..[[K]][[1]][[1]]   
      }
      
      message(" Step 1. Subsetting cases/columns: " , K, ".", k, " training set @ ", Sys.time())
      betas.train <- betas.mat[ , fold$train] # rows CpG # columns are patients! # see line 408
      
      if(subset.10K) {
        message("  Filtering top 10000 CpGs by row SD @ ", Sys.time())
        betas.train <- betas.train[order(apply(betas.train, 1, sd), decreasing = T)[1:10000], ] 
      }
      message("  Dimension of `betas.train` nrows: ", 
              nrow(betas.train), " ncols: ", ncol(betas.train))
      #        "\n  Variance filtering finished @ ", Sys.time()) # Duration @ single core ca. 1.25-1.5mins
      message(" \n  Check whether there is NA in train set : ", sum(is.na(betas.train) == T))
      
      betas.train <- t(betas.train)
      message("  Transposing `betas.train` finished @ ", Sys.time())
      
      message("  Clean up memory (garbage collector) @ ", Sys.time())
      gc()
      

      message(" Step 2. Subsetting `betas.mat` cases/columns: " , K, ".", k, " test/calibration set @ ", Sys.time())
      betas.test <- betas.mat[ , fold$test]
      
      if(subset.10K) {
        message("  Matching CpGs in test set to filtered train set @ ", Sys.time())
        betas.test <- betas.test[match(colnames(betas.train), rownames(betas.test)), ]
      }
      betas.test <- t(betas.test)
      message("  Transposing `betas.test` finished @ ", Sys.time())
      message("  Dimension of `betas.test`  nrows: ", 
              nrow(betas.test), " ncols: ", ncol(betas.test),
              "\n CpG matching finished @ ", Sys.time())
      
      if(!identical(colnames(betas.train), colnames(betas.test))) {
        warning(K,".",k, " Training/test set CpG column names do not match!")
      } else {
        message(K,".",k, " Training/test set CpG column names are consistent ✔")
      }
    
      folder.path <- file.path(getwd(), out.path)
      dir.create(folder.path, showWarnings = F, recursive = T)
      save(betas.train, 
           betas.test,
           fold,
           file = file.path(folder.path, paste(out.fname, K, k, "RData", sep = "."))
      )  
      message(K,".",k, " Results saved successfully @ ", Sys.time())
    }
  }
}


