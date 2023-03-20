#!/usr/bin/env R

# Author: Sean Maden
#
# Get differentially expressed genes (DEGs) among bulk sample experiment conditions.
#

source("deconvo_method-paper/code/11_bulk-expression-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
rse.markers <- get(load(rse.k2markers.filepath))

# compare mean expression across groups

# compare mean expression differences by groups
