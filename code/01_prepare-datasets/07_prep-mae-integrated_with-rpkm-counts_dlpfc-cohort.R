#!/usr/bin/env R

#
# Makes a multi assay experiment object containing the data for integrated analysis. 
# The new MAE object includes the following asssays:
# * snRNAseq
# * bulk RNAseq
# * RNAscope image processing outputs from HALO
#

libv <- c("here", "nlme", "ggplot2", "gridExtra", "dplyr", "ggforce", "MultiAssayExperiment", "SingleCellExperiment")
sapply(libv, library, character.only = TRUE)

#--------------------
# load prepped assays
#--------------------
# load snrnaseq data
load("~/GitHub/deconvo_method-paper/outputs/01_prepare-datasets/list-scef_markers-k2-k3-k4_ro1-dlpfc.rda")
# load bulk data
load("~/GitHub/deconvo_method-paper/outputs/01_prepare-datasets/rse-gene-filter.rda")
# load rnascope image data
load("~/GitHub/deconvo_method-paper/outputs/01_prepare-datasets/halo-outputs_updated.Rdata")
