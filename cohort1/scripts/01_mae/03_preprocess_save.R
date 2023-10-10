#!/usr/bin/env R

#
# Preprocess MultiAssayExperiment
#

filter.rnascope.confidence <- "Low"

#-----
# load
#-----

mae.in.path <- "./outputs/01_mae/mae_allsamples.rda"
mae <- get(load(mae.in.path))
cd <- colData(mae)

#-------------------------------------
# 2. filter on confidence annotations
#-------------------------------------
# filter on RNAscope confidence annotations

file.path <- "./data/09_quality/rnascope_quality_annotations.csv"
anno <- read.csv(file.path)
samples.vector.exclude <- c()

filter.star <- anno$Star=="Low"
samples.vector.exclude <- paste0(anno[filter.star,]$Sample, ";Star")

filter.circle <- anno$Circle=="Low"
samples.vector.exclude <- c(samples.vector.exclude, paste0(anno[filter.circle,]$Sample, ";Circle"))

filter.sce.img <- ""


sce.img$SAMPLE_ID

#-------------------------------------
# 3. filter on neuron proportion
#-------------------------------------

#-----
# save
#-----

mae.out.path <- "./outputs/01_mae/mae_analysis.rda"
mae <- get(load(mae.path))

