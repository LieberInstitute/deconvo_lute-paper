#!/usr/bin/env R

# Author: Sean Maden
#
# Get differentially expressed genes (DEGs) among bulk sample experiment conditions.
#

# load data
source("deconvo_method-paper/code/11_bulk-expression-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
rse.markers <- get(load(rse.k2markers.filepath))
halo <- get(load(halo.outputs.path))


