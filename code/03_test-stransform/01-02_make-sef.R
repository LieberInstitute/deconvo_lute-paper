#!/usr/bin/env R

#
# Make the filtered summarized experiment, non-DelayedArray.
# 

library(SingleCellExperiment)

#------
# paths
#------
proj.dpath <- "deconvo_method-paper"
save.dpath <- file.path(proj.dpath, "outputs/03_test-stransform")
sce.fpath <- "DLPFC_snRNAseq/processed-data/sce/sce_DLPFC.Rdata"
lz.fname <- "lz_s-rescale-k4_dlpfc-ro1.rda"

#-----
# load
#-----
# single cell experiment
sce <- get(load(sce.fpath))
# gene markers
lz <- get(load(file.path(save.dpath, lz.fname)))

#--------------
# remake as sef
#--------------