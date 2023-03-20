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

# correlate bulk and pseudobulk marker expression
correlation.matrix <- do.call(cbind, lapply(rse.experiment.groups, function(groupi){
  message(groupi)
  filter <- colData(rse.filter.markers)[[condition.variable]] == groupi
  rowMeans(assays(rse.filter.markers)[[assay.name]][,filter])
}))
colnames(correlation.matrix) <- rse.experiment.groups
# get correlation matrix
data.for.correction <- cbind(pseudobulk, correlation.matrix)
correlation.table <- cor(data.for.correction, method = correlation.method.markers)
# get correlation ggplot
ggplot.correlation.heatmap <- ggcorrplot(correlation.table, type = "lower", lab = T, 
                                         title = paste0(correlation.method.markers))
jpeg(correlation.heatmap.jpg.path, width = 10, height = 10, units = "in", res = 400)
print(ggplot.correlation.heatmap); dev.off()
