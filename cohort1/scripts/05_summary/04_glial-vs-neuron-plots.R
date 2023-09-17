#!/usr/bin/env R

# Author: Sean Maden
#
# Compare glial versus neuron cell size estimates.

source("deconvo_method-paper/code/07_cell-size-estimates/00_parameters.R")
sapply(libv, library, character.only = T)
sizes.table <- get(load(sizes.k2.table.path))

# correlation heatmap

# pairs plot

# glial vs. neuron