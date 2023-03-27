#!/usr/bin/env R

# Author: Sean Maden
#
# Compare cell sizes from snRNA-seq and image analysis.

source("deconvo_method-paper/code/07_cell-size-estimates/00_parameters.R")
sapply(libv, library, character.only = T)
sn.total <- get(load(sce.sizes.total.expression.path))
sn.gene <- get(load(sce.sizes.expressed.genes.path))
image.all <- get(load(image.cell.sizes.save.path))
# get harmonized matrix
# get glial/neuron
# save