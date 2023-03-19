#!/usr/bin/env R

# Author: Sean Maden
#
# Compare bulk, pseudobulk marker gene expression
#

# load data
source("deconvo_method-paper/code/11_bulk-expression-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
rse.markers <- get(load(rse.k2markers.filepath))
sce.markers <- get(load(sce.markers.list.path))[["k2"]]

# get pseudobulk from sce.markers
S <- c(3, 10)
P <- as.numeric(table(sce.markers[,"k2"]))

# compare expression -- counts

# compare expression -- means

# compare expression -- variances

# compare expression -- dispersion