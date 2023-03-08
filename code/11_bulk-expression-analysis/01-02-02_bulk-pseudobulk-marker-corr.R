#!/usr/bin/env R

# Author: Sean Maden
#
# Compare bulk, pseudobulk marker gene expression
#

libv <- c("SingleCellExperiment", "SummarizedExperiment", "glmGamPoi")
sapply(libv, library, character.only = T)

#----------
# load data
#----------
# get save directory path
save.path <- file.path("deconvo_method-paper", "outputs", 
                       "11_bulk-expression-analysis")

# load marker bulk expr
rsef.filename <- "rsef_k2-marker-expr_ro1-dlpfc.rda"
rsef <- get(load(file.path(save.path, rsef.filename)))

# load sce data
sce.markers.path <- file.path("deconvo_method-paper", "outputs", "09_manuscript", 
                              "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda")
lsce <- get(load(sce.markers.path))
sce <- lsce[["k2"]]
rm(lsce)

#-----------------------------
# compare expression -- counts
#-----------------------------

#-----------------------------
# compare expression -- means
#-----------------------------

#--------------------------------
# compare expression -- variances
#--------------------------------

#---------------------------------
# compare expression -- dispersion
#---------------------------------