#!/usr/bin/env R

# Author: Sean Maden
# 
# Validate sopt utility on new held-out bulk samples.
#

libv <- c("lute", "MultiAssayExperiment")
sapply(libv, library, character.only = T)

#-----
# load
#-----
# validation data
validation.sample.id <- c("Br6432", "Br6522", "Br8667")
# get all mae
mae.path <- file.path("deconvo_method-paper", "outputs", "01_prepare-datasets", "mae_with-rpkm_additional-data_final.rda")
mae <- get(load(mae.path))
# get sce marker data
sce <- mae[["sn1.rnaseq"]] # training marker expression
# subset mae
cd.mae <- colData(mae)
filter.validate.mae <- grepl(paste0(validation.sample.id,collapse = "|"),cd.mae$sample.id)
table(filter.validate.mae)
mae <- mae[,filter.validate.mae,]











