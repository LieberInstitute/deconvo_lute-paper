#!/usr/bin/env R

#
# Load and subset the SCE marker data for validation samples
#
#
libv <- c("SingleCellExperiment")
sapply(libv, library, character.only = T)
# load
base.path <- file.path("deconvo_method-paper", "outputs", "01_prepare-datasets")
# validation sce
sce.validate <- get(load(file.path(base.path, "sce-validate_k-2-3-4-markers_cohort1.rda")))
# training markers sce
list.sce.train.markers <- get(load(file.path(base.path, "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda")))
# subset validation sce into list
list.sce.validate.markers <- lapply(sce.train.markers, function(sce.iter){
  sce.validate[rownames(sce.validate) %in% rownames(sce.iter),]
})
names(list.sce.validate.markers) <- names(list.sce.train.markers)
# save list
save(list.sce.validate.markers, file = file.path(base.path, "list-sce-validate_markers-k-2-3-4_cohort1.rda"))