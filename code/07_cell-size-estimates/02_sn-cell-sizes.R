#!/usr/bin/env R

# Author: Sean Maden
#
# Get cell sizes from preprocessed snRNAseq data.

source("deconvo_method-paper/code/07_cell-size-estimates/00_parameters.R")
sapply(libv, library, character.only = T)
sce <- get(load(sce.path))
# get cell size estimates
list.total.expression.data <- sce_summary(sce, expression.summary.type = "total.expression")
list.expressed.genes.data <- sce_summary(sce, expression.summary.type = "expressed.genes")
# save
save(list.total.expression.data, file = sce.sizes.total.expression.path)
save(list.expressed.genes.data, file = sce.sizes.expressed.genes.path)