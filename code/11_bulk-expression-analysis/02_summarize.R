#!/usr/bin/env R

# Author: Sean Maden
#
# Summarize samples and genes.
#
source("deconvo_method-paper/code/11_bulk-expression-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
rse.filter <- get(load(rse.gene.filter.filepath))
cd <- colData(rse)
# summaries
table(cd$library_prep)
table(cd$library_type)
table(cd$library_type, cd$library_prep)
sample.variable <- paste0(cd$BrNum, "_", cd$location)
table(cd[,condition.variable])
