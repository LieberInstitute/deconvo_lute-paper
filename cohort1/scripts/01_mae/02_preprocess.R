#!/usr/bin/env R

#
# Preprocess MultiAssayExperiment
#

min.neuron.proportion <- 0.2
max.nucleus.area <- 78
filter.rnascope.confidence <- "Low"

#-----
# load
#-----

mae.in.path <- "./outputs/01_mae/mae_allsamples.rda"
mae <- get(load(mae.in.path))
cd <- colData(mae)

#--------------------------
# 1. filter on nucleus area
#--------------------------

assay.name <- "cell.sizes"

rnascope <- mae[[assay.name]]


dim(sce.img)
filter.sce <- assays(sce.img)[["Nucleus_Area"]] < max.nucleus.area
sce.img <- sce.img[,filter.sce]
dim(sce.img)

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

