#!/usr/bin/env R

# Author: Sean Maden
#
# Show the impact of including cell size adjustments on pseudobulking outcomes.

libv <- c("lute", "scater", "SingleCellExperiment", "SummarizedExperiment", "scater")
sapply(libv, library, character.only = T)

#----------
# load data 
#----------
read.dpath <- save.dpath <- file.path("deconvo_method-paper", "outputs", 
                                      "07_cell-size-estimates")

# cell size data
sce.csize.fname <- "df-cellsize_donor-region_sce.rda"
sce.csize <- get(load(file.path(read.dpath, sce.csize.fname)))

# sce marker data, k2
sef.fname <- "sef_mr-markers-k2-from-sce_dlpfc-ro1.rda"
sef.dpath <- file.path("deconvo_method-paper", "outputs", 
                        "05_marker-gene-annotations")
sef <- get(load(file.path(scef.dpath, scef.fname)))

#-------------------------
# simulations -- sce.csize
#-------------------------
# get cell sizes, s
dfs <- aggregate(sce.csize, by = list("celltype"), FUN = mean)

# get signature matrix, z
sce <- SingleCellExperiment(assays = list(counts = assays(scef)$co),
                            cd = colData(scef))
setf <- set_from_sce(scef, groupvar = "donor", method = "mean")



