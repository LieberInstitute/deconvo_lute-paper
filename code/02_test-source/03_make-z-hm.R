#!/usr/bin/env R

#
# Make some heatmaps from z datasets
# 

#------
# paths
#------
proj.dpath <- "deconvo_method-paper"
save.dpath <- file.path(proj.dpath, "outputs/02_test-source")
lz.fname <- "lz_mr_dlpfc-ro1.rda"

#-----
# load
#-----
lz <- get(load(file.path(save.dpath, lz.fname)))

#---------------
# table previews
#---------------


#---------
# heatmaps
#---------
pheatmap(lz$top.marker.data)




