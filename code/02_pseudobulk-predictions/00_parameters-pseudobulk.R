#!/usr/bin/env R

# Author: Sean Maden
#
# Main parameters, or dependency objects, for deconvolution framework trials.

# dependencies
libv <- c("here", "lute", "dplyr", "ggplot2", "gridExtra", "SingleCellExperiment", 
          "SummarizedExperiment", "scran")
sapply(libv, library, character.only = TRUE)

# save path
save.path <- here("deconvo_method-paper", "outputs", "02_pseudobulk-predictions")

# mrb sce path
sce.mrb.path <- here("deconvo_method-paper", "outputs", "01_prepare-datasets", "sce-mrb_dlpfc.rda")

# dlpfc markers path
sce.markers.list.path <- here("deconvo_method-paper", "outputs", "02_pseudobulk-predictions", "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda")