#!/usr/bin/env R

# Author: Sean Maden
#
# Run this script from: ./deconvo_method-paper/software/
#
# Makes multipanel plot of results from k2, k3, k4 experiments across cohorts, sn and bulk RNA-seq.
#
#

#-------------
# dependencies
#-------------

libv <- c("here", "nlme", "lute", "ggplot2", "gridExtra", "dplyr", "ggforce", "MultiAssayExperiment", "SingleCellExperiment")
sapply(libv, library, character.only = TRUE)

#-----------------
# assign variables
#-----------------
#

#-----
# save
#-----

# set new mae filename
new.mae.filename <- "mae_allsamples.rda"

mae.final.filepath <- file.path("outputs", "01_mae", new.mae.filename)

save(mae.final, file = mae.final.filepath)
