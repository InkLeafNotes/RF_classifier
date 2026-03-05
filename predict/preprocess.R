library(minfi)
source("MNPprocessIDAT_functions.R") ## from Github(https://github.com/mwsill/mnp_training)

##loading data
baseDir <- "."
targets <- read.metharray.sheet(baseDir)
idatdir <- "/media/user/data/Sanbo_Methylation/rawdata/"
targets$Basename <- paste(idatdir,"idat/",targets$Slide,"_",targets$Array,sep = "")

rgSet <- read.metharray.exp(targets=targets,force=T)
sampleNames(rgSet) <- targets$Sample_Name

## Illumina normalization
message("running normalization ...",Sys.time())
Mset <- MNPpreprocessIllumina(rgSet)

## batch effect
load("/media/user/data/GSE/GSE140686/preprocess/results/ba.coef.RData")

probes_Mset <- rownames(Mset)
probes_reference <- names(methy.coef$FFPE)
common_probes <- intersect(probes_Mset, probes_reference )
Mset <- Mset[common_probes, ]

batch <- rep("FFPE",length(Mset$Sample_Name))

methy <- getMeth(Mset)
unmethy <- getUnmeth(Mset)

# perform batch adjustment
methy.b <- log2(methy +1) + matrix(unlist(methy.coef[match(batch,names(methy.coef))]),ncol=length(batch))
unmethy.b <- log2(unmethy +1) + matrix(unlist(unmethy.coef[match(batch,names(unmethy.coef))]),ncol=length(batch))
methy.b[methy.b < 0] <- 0
unmethy.b[unmethy.b < 0] <- 0
methy.ba <- 2^methy.b
unmethy.ba <- 2^unmethy.b
# illumina-like beta values
betas.test <- methy.ba / (methy.ba +unmethy.ba +100)
betas.test <- as.data.frame(t(betas.test))
saveRDS(betas.test, "chordoma.ba.coef.rds")
