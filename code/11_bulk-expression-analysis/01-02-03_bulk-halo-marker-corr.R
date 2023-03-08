#!/usr/bin/env R

# Author: Sean Maden
#
# Get differentially expressed genes (DEGs) among bulk sample experiment conditions.
#

libv <- c("SingleCellExperiment", "SummarizedExperiment")
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

# load halo outputs

#--------------------------------------------
# correlation bulk and halo marker expression
#--------------------------------------------