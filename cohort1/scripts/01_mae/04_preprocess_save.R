#!/usr/bin/env R

#
# Preprocess MultiAssayExperiment
#

filter.rnascope.confidence <- "Low"

#-----
# load
#-----
mae.in.path <- "./outputs/01_mae/mae_allsamples_append.rda"
mae <- mae.all <- get(load(mae.in.path))
dim(colData(mae))

# rnascope confidence annotations
cd.id <- get(load("./outputs/01_mae/sample_qc_df.rda"))

#---------------------------------------------
# 1. filter on rnascope confidence annotations
#---------------------------------------------
sample.vector.keep <- cd.id[cd.id$remove.low==FALSE,]$sample.id

# filter sce.img
sce.img <- mae[["sce.img"]]
dim(sce.img)
filter.sce <- colData(sce.img)$Sample %in% sample.vector.keep
sce.img <- sce.img[,filter.sce]
dim(sce.img)
mae[["sce.img"]] <- sce.img

# filter cell.sizes
cell.sizes <- mae[["cell.sizes"]]
dim(cell.sizes)
filter.cell.sizes <- cell.sizes["sample_id",] %in% sample.vector.keep
cell.sizes <- cell.sizes[,filter.cell.sizes]
dim(cell.sizes)
mae[["cell.sizes"]] <- cell.sizes

#-----
# save
#-----

dim(colData(mae))
mae.out.path <- "./outputs/01_mae/mae_analysis_append.rda"
save(mae, file = mae.out.path)
