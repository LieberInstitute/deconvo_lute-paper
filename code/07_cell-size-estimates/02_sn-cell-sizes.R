#!/usr/bin/env R

# Author: Sean Maden
#
# Get cell sizes from preprocessed snRNAseq data.

source("deconvo_method-paper/code/11_bulk-expression-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
sce <- get(load(sce.path))