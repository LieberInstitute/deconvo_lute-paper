#!/usr/bin/env/R

# Author: Sean Maden
#
# Run K5 experiments on ABIS-seq data
#
#

libv <- c('tidyr')
sapply(libv, library, character.only = TRUE)

#-----
# load
#-----

load("./env/02_abisseq/01_abisseq_script.RData")
load("./env/02_abisseq/02_proportions_s13_script.RData")
source("./source/lute_experiment.R")
source("./source/cell_mappings_helpers.R")

#------------------
# map k5 cell types
#------------------

vectorCellTypeMap <- c("T", "B", "Dendrocyte", "Plasma", "Monocyte", "Neutrophil", "NK")
vectorCellTypeStart <- colnames(zref)

mappingsTable <- cellLabelMappings(
  vectorCellTypeMap, vectorCellTypeStart, returnType = "list")[[1]]

zrefMapped <- cellLabelMappings(
  vectorCellTypeMap, vectorCellTypeStart, zref, 
  returnType = "df", summaryOperation = "mean")

trueProportionsStart <- fc.proportions
sampleNames <- fc.proportions[,1]
trueProportionsStart <- trueProportionsStart[
  ,colnames(trueProportionsStart) %in% colnames(zref)]
rownames(trueProportionsStart) <- sampleNames
trueProportionsStart <- trueProportionsStart[,c(2:ncol(trueProportionsStart))]
trueproportionsMapped <- cellLabelMappings(
  vectorCellTypeMap, colnames(trueProportionsStart), 
  trueProportionsStart, returnType = "df", summaryOperation = "sum")
# remove nas
trueproportionsMapped <- trueproportionsMapped[
  ,!colnames(trueproportionsMapped)=="NA"]
# get fractions <1
trueproportionsMapped <- t(apply(
  trueproportionsMapped, 1, function(ri){ri/sum(ri)}))
trueproportionsMapped <- as.data.frame(trueproportionsMapped)
rownames(trueproportionsMapped) <- sampleNames

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

bulkExpression <- as.matrix(assays(se)[["tpm"]])
bulkExpression <- bulkExpression[rownames(bulkExpression) %in% rownames(zref),]
referenceExpression <- zrefMapped
referenceExpression <- referenceExpression[
  rownames(referenceExpression) %in% rownames(bulkExpression),]
trueCellTypeProportions <- trueproportionsMapped

# get experiment results
experimentList <- newExperimentList(
  referenceExpression=referenceExpression,
  trueCellTypeProportions=trueCellTypeProportions,
  bulkExpression=bulkExpression,
  typemarkerAlgorithmName=NULL
)

experimentResults <- evaluateExperiment(experimentList, TRUE)
  
  
  
  
  singleCellExperiment=NULL,
  referenceExpression=zrefMapped,
  bulkExpression=tpm,
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

z = as.matrix(zref), 
y = as.matrix(assays(se)[["tpm"]]), 
assay.name = 'tpm',
typemarker.algorithm = NULL

#-----------
# save image
#-----------
save(file="env/03_k5/01_run_script.RData")
