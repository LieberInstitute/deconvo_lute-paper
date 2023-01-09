#!/usr/bin/env R

# Author: Sean Maden
#
# Show the impact of including cell size adjustments on pseudobulking outcomes.

libv <- c("lute", "SingleCellExperiment", "SummarizedExperiment")
sapply(libv, library, character.only = T)

#----------
# load data 
#----------
read.dpath <- save.dpath <- file.path("deconvo_method-paper", "outputs", 
                                      "07_cell-size-estimates")

# cell size data
sce.csize.fname <- "df-cellsize_donor-region_sce.rda"
sce.csize <- get(load(file.path(read.dpath, sce.csize.fname)))

# se marker data, k2
sef.fname <- "sef_mr-markers_k2_20-per-k_dlpfc-ro1.rda"
sef.dpath <- file.path("deconvo_method-paper", "outputs", 
                        "05_marker-gene-annotations")
sef <- get(load(file.path(sef.dpath, sef.fname)))

#-------------------------
# simulations -- sce.csize
#-------------------------
# get cell sizes, s
dfs <- aggregate(sce.csize, by = list("celltype"), FUN = mean)
# get signature matrix, z
sef[["donor"]] <- sef[["BrNum"]]
setf <- set_from_sce(sef, groupvar = "donor", method = "mean",
                     assayname = "logcounts")
