#!/usr/bin/env R

# Author: Sean Maden
#
#
source("deconvo_method-paper/code/11_bulk-expression-analysis/00_parameters.R")
sapply(libv, library, character.only = T)

# load data
rse.background <- get(load(rse.bulk.filepath))
rse.markers <- get(load(rse.k2markers.filepath))

# prepare comparison params
# set up data summaries
# params

# rse data
cd <- colData(rse)
# get new cd variables
cd[,batch.variable] <- paste0(cd$BrNum,"_",cd$location)
cd[,condition.variable] <- paste0(cd$library_prep,"_",cd$library_type)

# get expression scales to analyze
lcompare <- lapply(assays, function(assay.name){
  expression.background <- assays(rse)[[assay.name]]
  expression.markers <- assays(rsef)[[assay.name]]
  get_comparison_data(expr.bg = counts.bg, expr.marker = counts.marker,
                      plot.fname = plot.fname, type.vector = type.vector,
                      cd = cd, save.path = save.path)
})
