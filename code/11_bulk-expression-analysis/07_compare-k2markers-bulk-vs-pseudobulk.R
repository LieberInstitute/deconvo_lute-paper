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
rse.experiment.groups <- unique(rse.filter.markers[[condition.variable]])

# compare expression
# means
# variances
# dispersion

# save plots
jpeg(correlation.heatmap.jpg.path, width = 10, height = 10, units = "in", res = 400)
print(ggplot.correlation.heatmap); dev.off()
