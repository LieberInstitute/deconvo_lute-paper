#!/usr/bin/env R

# Author: Sean Maden
#
# Get differentially expressed genes (DEGs) among bulk sample experiment conditions.
#

source("deconvo_method-paper/code/11_bulk-expression-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
rse.markers <- get(load(rse.k2markers.filepath))

# perform differential expression analysis
se.markers <- SummarizedExperiment(rse.markers)
cd <- colData(rse.markers)[,condition.variable,drop=F]
names(assays(se.markers)) <- "counts"
des.markers <- DESeqDataSet(se.markers, design = cd)
count.data <- as.matrix(assays(rse.markers)[[1]])
des.markers <- DESeqDataSetFromMatrix(countData = count.data, 
                                      colData = cd, design = design.matrix)
# pairwise marker degs
experiment.group.variable <- "experiment.group"
experiment.group.vector <- colData(rsef[,experiment.group.variable])
unique.groups <- unique(experiment.group.vector)
combinations <- combn(unique.groups, 2)
dds <- DESeq(dds)
res <- results(dds)
