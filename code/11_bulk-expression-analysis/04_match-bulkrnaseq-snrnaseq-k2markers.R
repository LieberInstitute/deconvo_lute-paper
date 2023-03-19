#!/usr/bin/env R

# Author: Sean Maden
#
#

source("deconvo_method-paper/code/11_bulk-expression-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
rse.filter <- get(load(rse.gene.filter.filepath))
sce <- get(load(sce.markers.list.path))[["k2"]]

# get bulk marker expression
marker.genes.vector <- rownames(sce)
# subset rse on markers
bulk.gene.names <- rowData(rse.filter)$Symbol
gene.marker.intersect <- intersect(bulk.gene.names, marker.genes.vector)
message("Found ", length(gene.marker.intersect), " overlapping markers.")
rse.filter.markers <- which(bulk.gene.names %in% gene.marker.intersect)
rse.filter.markers <- rse.filter[rse.filter.markers,]
# gene symbols as rownames
rownames(rse.filter.markers) <- rowData(rse.filter.markers)$Symbol
dim(rse.filter.markers)
# save
save(rse.filter.markers, file = rse.k2markers.filepath)

# get pseudobulk marker expression
pseudobulk <- get_pseudobulk(sce, sce.assay.name, S)
# match markers
order.pseudobulk <- order(match(rownames(pseudobulk), 
                                rownames(rse.filter.markers)))
pseudobulk.ordered <- as.matrix(pseudobulk[order.pseudobulk,], ncol = 1)
if(!identical(rownames(pseudobulk.ordered), rownames(rse.filter.markers))){
  stop("Error, couldn't match markers in pseudobulk and bulk/rse")} else{
    pseudobulk <- pseudobulk.ordered
  }
colnames(pseudobulk) <- "pseudobulk"
# save
save(pseudobulk, file = pseudobulk.path)

# correlate bulk and pseudobulk marker expression
# get correlation matrix
rse.marker.expression <- assays(rse.filter.markers)[[assay.name]]
rse.experiment.groups <- unique(rse.filter.markers[[condition.variable]])
correlation.matrix <- do.call(cbind, lapply(rse.experiment.groups, function(groupi){
  message(groupi)
  filter <- colData(rse.filter.markers)[[condition.variable]] == groupi
  rowMeans(assays(rse.filter.markers)[[assay.name]][,filter])
}))
colnames(correlation.matrix) <- rse.experiment.groups
# get correlation matrix
correlation.matrix <- cbind(pseudobulk, correlation.matrix)
plot.data.frame <- cor(correlation.matrix, method = correlation.method.markers)
# get correlation ggplot
ggplot.correlation.heatmap <- ggcorrplot(plot.data.frame, type = "lower", lab = T, 
                                         title = paste0(correlation.method.markers))
jpeg(correlation.heatmap.jpg.path, width = 10, height = 10, units = "in", res = 400)
print(ggplot.correlation.heatmap); dev.off()
