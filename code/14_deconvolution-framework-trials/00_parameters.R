#!/usr/bin/env R

# Author: Sean Maden
#
# Main parameters, or dependency objects, for deconvolution framework trials.

# dependencies
libv <- c("here", "lute", "dplyr", "ggplot2", "gridExtra", "SingleCellExperiment", 
          "SummarizedExperiment")
sapply(libv, library, character.only = TRUE)

# helper functions

# load cell sizes
cell.sizes.k2.path <- here("deconvo_method-paper", "outputs", "07_cell-size-estimates")
cell.sizes.k2.path <- here(cell.sizes.k2.path, "cell-sizes-k2-table.rda")

#------------------
# script parameters
#------------------
# 01, independent pseudobulk
sce.mrb.name <- "list-scef_markers-k2-k3-k4_mrb-dlpfc"
sce.mrb.path <- here("deconvo_method-paper", "outputs", "09_manuscript", sce.mrb.name)

# 02, within-samples tests

# 03, across-samples tests

# 04