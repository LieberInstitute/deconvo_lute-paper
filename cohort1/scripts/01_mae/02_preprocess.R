#!/usr/bin/env R

#
# Preprocess MultiAssayExperiment
#

min.neuron.proportion <- 0.2
max.nucleus.area <- 78

#-----
# load
#-----

mae.in.path <- "./outputs/01_mae/mae_allsamples.rda"
mae <- get(load(mae.in.path))

#----------------------------------------------
# 1. filter sample IDs on cell type proportions
#----------------------------------------------

#------------------------------------------
# 2. filter RNAscope/HALO cells on max area
#------------------------------------------

#-----
# save
#-----

mae.out.path <- "./outputs/01_mae/mae_analysis.rda"
mae <- get(load(mae.path))

