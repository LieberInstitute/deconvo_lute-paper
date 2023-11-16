#!/usr/bin/env R

# Author: Sean Maden
#
# Read in the experiment elements.
#


processZ <- function(path){
  # processZ
  #
  # Get reference from file.
  
  # load
  zref <- 
    read.csv(
      "./data/monaco_et_al_2019/manuscript/abisseq_rnaseq_cell_types_references_k17.csv")
  
  # process
  zref <- zref[!zref[,1]=="12:00 AM",]
  rownames(zref) <- zref[,1]
  zref <- zref[,c(2:ncol(zref))]
  
  return(zref)
}

processP <- function(path){
  # processP
  # Get proportions from file.
  #
  
  p <- read.csv(path)
  
  # process
  rownames(p) <- p[,1]
  p <- p[,seq(2,ncol(p))]
  return(p)
}

processY <- function(path, assayName = "tpm"){
  # processY
  #
  # Get bulk expression from file.
  
  # load
  yTable <- read.table(path)
  
  
  
  #---------
  
  
  # process
  
  
  #---------
  assayList <- list(yTable)
  names(assayList) <- assayName
  ySe <- SummarizedExperiment(assays = assayList)
  returnList <- list(
    yTable=yTable,
    ySe=ySe
  )
  
  # get coldata
  # get s13 phenotype info from column labels
  phenoDataS13 <- colnames(yTable)
  table(gsub(".*_", "", phenoDataS13))
  phenoDataS13 <- data.frame(sample.id = phenoDataS13,
                             source.id = gsub("_.*", "", phenoDataS13),
                             sample.type = gsub(".*_", "", phenoDataS13))
  phenoDataS13$tissue.type <- 
    ifelse(phenoDataS13$sample.type=="PBMC", "PBMC", "immune_cell")
  phenoDataS13$tissue.type.detail <- 
    ifelse(phenoDataS13$sample.type=="PBMC", "PBMC", phenoDataS13$sample.type)
  # append summary statistics to pheno data
  yTable <- as.matrix(yTable)
  phenoDataS13$library.size <- colSums(yTable)
  phenoDataS13$mean.expression <- colMeans(yTable)
  phenoDataS13$median.expression <- colMedians(yTable)
  phenoDataS13$sd.expression <- colSds(yTable)
  phenoDataS13$num.na.expression <- colAnyNAs(yTable)
  phenoDataS13$num.zero.expression <- unlist(apply(yTable, 2, function(ci){length(ci[ci==0])}))
  rownames(phenoDataS13) <- phenoDataS13[,1]
  # append pheno data to colData for SummarizedExperiment
  colData(newSummarizedExperiment) <- DataFrame(phenoDataS13)
  
  
  return(returnList)
}


getExperimentData <- function(){
  # getExperimentData
  #
  # Get list of experiment data from paths.
  
  # get Y
  yPath <- "./data/monaco_et_al_2019/geo/GSE107011_Processed_data_TPM.txt"
  yList <- processY(yPath)
  
  # get Z
  zPath <- "./data/monaco_et_al_2019/manuscript/abisseq_rnaseq_cell_types_references_k17.csv"
  zref <- processZ(zPath)
  
  # get P
  pPath <- "./data/monaco_et_al_2019/manuscript/flow_cytometry_true_cell_type_proportions_s13.csv"
  ptrue <- processP(pPath)
  
  returnList <- list(
    y.table = yList[["yTable"]],
    y.se = yList[["ySe"]],
    z = zref,
    p.true = ptrue
  )
  return(returnList)
}