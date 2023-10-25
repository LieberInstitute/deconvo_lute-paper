#!/usr/bin/env/R

# Author: Sean Maden
#
# Run K5 experiments on ABIS-seq data
#
#

#-----
# load
#-----
load("./env/02_abisseq/01_abisseq_script.RData")
source("./source/lute_experiments.R")

#------------------
# map k5 cell types
#------------------
cellTypes <- c("T", "B", "Dendrocyte", "Plasma", "Monocyte", "Neutrophil", "NK")
cellTypesRemove <- c("")
cellSizes <- c()

#---------------
# run experiment
#---------------
# get experiment results
experimentList <- newExperimentList(
  singleCellExperiment=singleCellExperimentDarmanis,
  referenceExpression=NULL,
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
