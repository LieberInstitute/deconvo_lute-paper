#!/usr/bin/env R

# Author: Sean Maden
#
#

libv <- c("SingleCellExperiment", "SummarizedExperiment")
sapply(libv, library, character.only = T)

#----------
# load data
#----------
# load bulk data
rse.filename <- "rse_gene.Rdata"
load.path <- file.path("Human_DLPFC_Deconvolution",
                       "processed-data",
                       "01_SPEAQeasy",
                       "round2_v40_2022-07-06",
                       "rse")
rse <- get(load(file.path(path, rse.filename)))

# get save directory path
save.path <- file.path("deconvo_method-paper", "outputs", 
                       "11_bulk-expression-analysis")

# load marker data
sce.markers.path <- file.path("deconvo_method-paper", "outputs", "09_manuscript", 
                              "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda")
lsce <- get(load(sce.markers.path))

#----------------------------
# marker gene bulk expression
#----------------------------
