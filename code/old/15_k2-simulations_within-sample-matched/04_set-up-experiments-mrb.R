#!/usr/bin/env R

# Author: Sean Maden
#

source("deconvo_method-paper/code/15_sample-wise-signature-matrix-simulations/00_parameters.R")
sapply(libv, library, character.only = T)
sce.mrb <- get(load(sce.mrb.path))

# get lognorm counts
sce.mrb <- logNormCounts(sce.mrb, assay.type = "counts")

# get k2 labels and filter
mrb.cell.type.vector <- sce.mrb[["cellType"]]
mrb.is.neuron <- grepl("Excit|Inhib", mrb.cell.type.vector)
mrb.is.glial <- grepl("Oligo|OPC|Micro|Astro", mrb.cell.type.vector)
colData(sce.mrb)[,"k2"] <- ifelse(mrb.is.neuron, "neuron", ifelse(mrb.is.glial, "glial", "other"))
filter.cell.type <- sce.mrb[["k2"]] %in% c("neuron", "glial")
sce.mrb <- sce.mrb[,filter.cell.type]

# get markers by batch
mrb.markers.by.batch <- markers_by_batch(sce = sce.mrb, "donor", "k2", "logcounts", 20)

# plot overlaps
list.markers <- lapply(mrb.markers.by.batch, function(markers){markers$gene})
upset(fromList(list.markers))
