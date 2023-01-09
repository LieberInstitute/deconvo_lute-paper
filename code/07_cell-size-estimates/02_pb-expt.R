#!/usr/bin/env R

# Author: Sean Maden
#
# Show the impact of including cell size adjustments on pseudobulking outcomes.

libv <- c("lute")
sapply(libv, library, character.only = T)

#----------
# load data 
#----------
read.dpath <- save.dpath <- file.path("deconvo_method-paper", "outputs", 
                                      "07_cell-size-estimates")

# cell size data
sce.csize.fname <- "df-cellsize_donor-region_sce.rda"
sce.csize <- get(load(file.path(read.dpath, sce.csize.fname)))

#-------------------------
# simulations -- sce.csize
#-------------------------
dfs <- aggregate(sce.csize, by = list("celltype"), FUN = mean)



