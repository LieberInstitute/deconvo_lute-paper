#!/usr/bin/env R

# Author: Sean Maden
#
# Read in the experiment elements.
#



getExperimentData <- function(){
  # get Y
  geoTpmTable <- "./data/monaco_et_al_2019/geo/GSE107011_Processed_data_TPM.txt"
  yTable <- read.table(geoTpmTable)
  assayName <- "tpm"
  assayList <- list(item1=yTable)
  names(assayList) <- assayName
  ySe <- SummarizedExperiment(assays = assayList)
  
  # get Z
  zref <- 
    read.csv(
      "./data/monaco_et_al_2019/manuscript/abisseq_rnaseq_cell_types_references_k17.csv")
  rownames(zref) <- zref[,1]
  zref <- zref[,seq(2:ncol(zref))]
  
  # get P
  ptrue <- 
    read.csv(
      "./data/monaco_et_al_2019/manuscript/flow_cytometry_true_cell_type_proportions_s13.csv")
  
  returnList <- list(
    y.table = yTable,
    y.se = ySe,
    z = zref,
    p.true = ptrue
    
  )
  return(returnList)
}