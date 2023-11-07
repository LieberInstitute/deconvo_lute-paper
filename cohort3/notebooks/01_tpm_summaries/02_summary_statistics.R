#!/usr/bin/env R

# Author: Sean Maden
#
# Appends summary statistics to the SummarizedExperiment object.
#
#
#
#

libv <- c("SummarizedExperiment", "biomaRt")
sapply(libv, library, character.only = TRUE)

# save new SummarizedExperiment
summarizedExperimentPath <- "outputs/01_tpm_summaries/se_tpm_s13.rda"
summarizedExperiment <- get(load(summarizedExperimentPath))
phenoDataS13 <- colData(summarizedExperiment)

# append summary statistics to pheno data
for(assayType in names(assays(summarizedExperiment))){
  
  expressionTable <- as.matrix(assays(summarizedExperiment)[[assayType]])
  phenoDataS13[,ncol(phenoDataS13)+1] <- colSums(expressionTable)
  colnames(phenoDataS13)[ncol(phenoDataS13)] <- 
    paste0("library.size.", assayType)
  
  phenoDataS13[,ncol(phenoDataS13)+1] <- colMeans(expressionTable)
  colnames(phenoDataS13)[ncol(phenoDataS13)] <- 
    paste0("mean.expression.", assayType)
  
  phenoDataS13[,ncol(phenoDataS13)+1] <- colMedians(expressionTable)
  colnames(phenoDataS13)[ncol(phenoDataS13)] <- 
    paste0("median.expression.", assayType)
  
  phenoDataS13[,ncol(phenoDataS13)+1] <- colSds(expressionTable)
  colnames(phenoDataS13)[ncol(phenoDataS13)] <- 
    paste0("sd.expression.", assayType)
  
  phenoDataS13[,ncol(phenoDataS13)+1] <- colAnyNAs(expressionTable)
  colnames(phenoDataS13)[ncol(phenoDataS13)] <- 
    paste0("num.na.expression.", assayType)
  
  phenoDataS13[,ncol(phenoDataS13)+1] <- 
    unlist(
      apply(
        tpm, 2, function(ci){length(ci[ci==0])}))
  colnames(phenoDataS13)[ncol(phenoDataS13)] <- 
    paste0("num.zero.expression.", assayType)
  
  
  
  
}

#-----
# save
#-----

# save SummarizedExperiment
save(
  newSummarizedExperiment, 
  file = "outputs/01_tpm_summaries/se_tpm_s13_colData_summaries.rda")

# save env
save.image("./env/01_tpm_summaries/02_summary_statistics_script.RData")
