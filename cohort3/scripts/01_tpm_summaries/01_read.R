#!/usr/bin/env R

# Author: Sean Maden
#
# Read in TPM expression data.

libv <- c("SummarizedExperiment")
sapply(libv, library, character.only = TRUE)

#-----
# load
#-----
tpm.path <- "./data/GSE107011_Processed_data_TPM/GSE107011_Processed_data_TPM.txt"
tpm <- read.table(tpm.path)

#------------------
# get pheno/coldata
#------------------
cd <- colnames(tpm)
table(gsub(".*_", "", cd))
cd <- data.frame(sample.id = cd, 
                 source.id = gsub("_.*", "", cd),
                 sample.type = gsub(".*_", "", cd))
cd$tissue.type <- ifelse(cd$sample.type=="PBMC", "PBMC", "immune_cell")
cd$tissue.type.detail <- ifelse(cd$sample.type=="PBMC", "PBMC", cd$sample.type)
# append summary statistics
tpm <- as.matrix(tpm)
cd$library.size <- colSums(tpm)
cd$mean.expression <- colMeans(tpm)
cd$median.expression <- colMedians(tpm)
cd$sd.expression <- colSds(tpm)
cd$num.na.expression <- colAnyNAs(tpm)
cd$num.zero.expression <- unlist(apply(tpm, 2, function(ci){length(ci[ci==0])}))
rownames(cd) <- cd[,1]
# format summarized experiment
se <- SummarizedExperiment(assays = list(tpm = tpm))
colData(se) <- DataFrame(cd)

#-----
# save
#-----
save.image("./env/01_tpm_summaries/01_read_script.RData")