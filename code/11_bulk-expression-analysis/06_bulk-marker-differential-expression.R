#!/usr/bin/env R

# Author: Sean Maden
#
# Get differentially expressed genes (DEGs) among bulk sample experiment conditions.
#

source("deconvo_method-paper/code/11_bulk-expression-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
rse.markers <- get(load(rse.k2markers.filepath))

# pairwise marker degs
experiment.group.variable <- "experiment.group"
experiment.group.vector <- colData(rsef[,experiment.group.variable])
unique.groups <- unique(experiment.group.vector)
combinations <- combn(unique.groups, 2)
dds <- DESeq(dds)
res <- results(dds)
