#!/usr/bin/env/R

# Author: Sean Maden
#
# Run K5 experiments on ABIS-seq data
#
#

library(tidyr)

#-----
# load
#-----

load("./env/02_abisseq/01_abisseq_script.RData")
source("./source/lute_experiments.R")

#------------------
# map k5 cell types
#------------------

vectorCellTypeMap <- c("T", "B", "Dendrocyte", "Plasma", "Monocyte", "Neutrophil", "NK")
vectorCellTypeStart <- colnames(zref)

mappingsTable <- cellLabelMappings(
  vectorCellTypeMap, vectorCellTypeStart, returnType = "list")[[1]]

zrefMapped <- cellLabelMappings(
  vectorCellTypeMap, vectorCellTypeStart, zref, returnType = "df")

# tests
#
#
#
# test label Monocyte
mappedLabel <- "Monocyte"
# test gene FXYD6, label Monocyte
geneName <- "FXYD6"
zrefMapped[geneName,mappedLabel]==rowMeans(
  zref[,which(colnames(zref) %in% mappingsTable[mappingsTable[,2]==mappedLabel,1])]
)[geneName]
# test gene NRG1, label Monocyte
geneName <- "NRG1"
zrefMapped[geneName,mappedLabel]==rowMeans(
  zref[,which(colnames(zref) %in% mappingsTable[mappingsTable[,2]==mappedLabel,1])]
)[geneName]

# test label Plasma
mappedLabel <- "Plasma"
# test gene FXYD6, label Plasma
geneName <- "FXYD6"
zrefMapped[geneName,mappedLabel]==
  zref[,which(colnames(zref) %in% 
                mappingsTable[mappingsTable[,2]==mappedLabel,1])][
                  which(rownames(zref)==geneName)]
# test gene NRG1, label Plasma
geneName <- "NRG1"
zrefMapped[geneName,mappedLabel]==
  zref[,which(colnames(zref) %in% 
                mappingsTable[mappingsTable[,2]==mappedLabel,1])][
                  which(rownames(zref)==geneName)]








#---------------
# run experiment
#---------------
# get experiment results
experimentList <- newExperimentList(
  singleCellExperiment=NULL,
  referenceExpression=zref.new,
  bulkExpression=bulkExpressionDarmanis,
  cellScaleFactors=cellScaleFactors,
  experimentType=experimentType,
  trueCellTypeProportions=
    trueCellTypeProportions,
  trueCellTypeProportionsSource=
    trueCellTypeProportionsSource,
  assayName="counts",
  cellTypeVariable="celltype",
  deconvolutionAlgorithmName="nnls",
  typemarkerAlgorithmName=NULL)
experimentResults <- evaluateExperiment(experimentList, TRUE)

#-----------
# save image
#-----------
save(file="env/03_k5/01_run_script.RData")
