#!/usr/bin/env R

#
# Preprocess MultiAssayExperiment
#

filter.rnascope.confidence <- "Low"

#-----
# load
#-----

mae.in.path <- "./outputs/01_mae/mae_allsamples_append.rda"
mae <- get(load(mae.in.path))
cd <- colData(mae)

#---------------------------------------------
# 1. filter on rnascope confidence annotations
#---------------------------------------------
# sample ids to remove
sce.img <- mae[["sce.img"]]
filter.remove <- sce.img$Confidence=="Low"
sample.id.remove <- unique(sce.img[,filter.remove]$Sample)
filter.mae <- !colData(mae)$sample.id %in% sample.id.remove
mae <- mae[,filter.mae,]

#-----
# save
#-----
mae.out.path <- "./outputs/01_mae/mae_analysis_append.rda"
mae <- get(load(mae.out.path))
