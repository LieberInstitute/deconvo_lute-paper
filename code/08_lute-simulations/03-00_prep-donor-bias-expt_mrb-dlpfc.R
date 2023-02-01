#!/usr/bin/env R

# Author: Sean Maden
#
# Prepare datasets for donor bias experiments using the multi-region brain 
# DLPFC dataset.
#
#

libv <- c("scuttle", "glmGamPoi", "sva",
          "SingleCellExperiment", "SummarizedExperiment")
sapply(libv, library, character.only = TRUE)


#----------
# load data
#----------
# get save dpath
code.dname <- "08_lute-simulations"
proj.dname <- "deconvo_method-paper"
save.dpath <- file.path(proj.dname, "outputs", code.dname)
# get marker data
dfm.fname <- "markers-k2_db-mr2_sce-dlpfc-mrb.rda"
dfm <- get(load(file.path(save.dpath, dfm.fname)))
# get sce
sce.fname <- "sce-mrb_dlpfc.rda"
sce.fpath <- file.path(save.dpath, sce.fname)
sce <- get(load(sce.fpath))

#-------------------
# get top 20 markers
#-------------------

#---------------------
# perform downsampling
#---------------------

# save

#---------------
# perform combat
#---------------

# save

#-----------------------------
# get neg. binom. coefficients
#-----------------------------
# summarize by donor

# fit models

# get model coeffs

# save