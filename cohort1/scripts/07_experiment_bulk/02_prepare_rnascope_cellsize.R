#!/usr/bin/env R

#
# Gets the RNAscope cell size estimate from MAE data.
#
#

libv <- c("ggplot2", "dplyr", "MultiAssayExperiment")
sapply(libv, library, character.only = T)

# load mae (SEE CODE 01 OUTPUTS)
new.mae.filename <- "mae_with-rpkm_additional-data_final.rda"
mae.final.filepath <- file.path("outputs", "01_prepare-datasets", new.mae.filename)
mae <- get(load(mae.final.filepath))

# Get an RNAscope cell size estimate
df.ct.info <- mae[["df.cellstat.rnascope"]]
neuron.size <- median(as.numeric(df.ct.info["cell_size",grepl("neuron", colnames(df.ct.info))]))
glial.size <- median(as.numeric(df.ct.info["cell_size",grepl("glial", colnames(df.ct.info))]))