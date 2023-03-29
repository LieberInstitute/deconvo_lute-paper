#!/usr/bin/env R

# Author: Sean Maden
#
# This script performs the deconvolution trials within samples, using matched
# snRNAseq, bulk RNAseq, and image outputs.

source("deconvo_method-paper/code/14_deconvolution-framework-trials/00_parameters.R")
sapply(libv, library, character.only = T)
rse <- get(load(rse.gene.filter.filepath))
image.table <- get(load(halo.output.path))
sce <- get(load(sce.markers.list.path))[["k2"]]

# get matched datasets
# image
image.sample.id.vector <- image.table[,"BrNum"]
image.sample.id.vector <- paste0(image.sample.id.vector, "_",
                                 toupper(image.table[,"Position"]))
# sn
sce.sample.id.vector <- sce[["BrNum"]]
sce.sample.id.vector <- paste0(sce.sample.id.vector, "_",
                               toupper(sce[["Position"]]))
# bulk
rse.sample.id.vector <- rse[["BrNum"]]
rse.sample.id.vector <- paste0(rse.sample.id.vector, "_", 
                               toupper(rse[["location"]]))









