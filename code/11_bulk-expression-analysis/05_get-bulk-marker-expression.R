#!/usr/bin/env R

# Author: Sean Maden
#
#

# load data
source("deconvo_method-paper/code/11_bulk-expression-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
rse.background <- get(load(rse.gene.filter.filepath))
rse.markers <- get(load(rse.gene.filter.filepath))

# prepare comparison params
# set up data summaries
# params

# get expression scales to analyze
lcompare <- lapply(assays, function(assay.name){
  expression.background <- assays(rse.background)[[assay.name]]
  expression.markers <- assays(rse.markers)[[assay.name]]
  get_comparison_data(expr.bg = counts.bg, expr.marker = counts.marker,
                      plot.fname = plot.fname, type.vector = type.vector,
                      cd = cd, save.path = save.path)
})


# qc at markers versus background -- jitter/boxplots
# mean expression
# expression variance
# correlation bulk and halo marker expression
# marker degs among bulk samples