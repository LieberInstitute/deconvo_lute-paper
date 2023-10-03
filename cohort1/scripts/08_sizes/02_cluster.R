#!/usr/bin/env R

# Author: Sean Maden
#
# Get cell type proportions
#

libv <- c("dplyr", "ggplot2")
sapply(libv, library, character.only = TRUE)
folder.name <- "08_sizes"

#------
# load
#------
load(file.path("./env/",folder.name,"/01_sizes_script.RData"))

#-----------------------------------
# pca cluster sizes, color by region
#-----------------------------------
dfp.tall$neuron
prcomp()

#--------
# save
#--------
# env
save.image(file.path("./env/",folder.name,"/02_cluster_script.RData"))
