#!/usr/bin/env R

# Author: Sean Maden
#

source("deconvo_method-paper/code/15_sample-wise-signature-matrix-simulations/00_parameters-script-set-15.R")
sapply(libv, library, character.only = T)
sce <- get(load(sce.path))

# get lognorm counts
sce <- logNormCounts(sce, assay.type = "counts")

# get k2 labels and filter
cell.type.vector <- sce[["cellType_broad_hc"]]
k2.cell.type.vector <- ifelse(cell.type.vector %in% c("Excit", "Inhib"), "neuron",
                              ifelse(cell.type.vector %in% c("Oligo", "OPC", "Astro", "Micro"), 
                                     "glial", "other"))
sce[["k2"]] <- k2.cell.type.vector
sce <- sce[,!sce[["k2"]]=="other"]
save(sce, file = sce.prepared.path)

# get markers by batch
# markers.by.batch <- markers_by_batch(sce = sce.mrb, "donor", "k2", "logcounts", 20)

# visualize overlaps
# list.markers <- lapply(markers.by.batch, function(markers){markers$gene})
# upset(fromList(list.markers))
