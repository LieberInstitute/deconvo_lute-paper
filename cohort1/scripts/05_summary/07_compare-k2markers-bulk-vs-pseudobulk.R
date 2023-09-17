#!/usr/bin/env R

# Author: Sean Maden
#
# Compare bulk, pseudobulk marker gene expression
#

# load data
source("deconvo_method-paper/code/11_bulk-expression-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
rse.markers <- get(load(rse.k2markers.filepath))
pseudobulk <- get(load(pseudobulk.path))
# get bulk experiment groups
rse.marker.expression <- assays(rse.filter.markers)[[assay.name]]
rse.experiment.groups <- rse.filter.markers[[condition.variable]] %>% unique()

# compare expression
# means
# variances
# dispersion

# save plots
jpeg(, width = , height = , units = "in", res = 400)
print(); dev.off()
