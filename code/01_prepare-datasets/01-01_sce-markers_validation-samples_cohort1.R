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

# append k marker labels
list.sce.validate.markers <- lapply(list.sce.validate.markers, function(sce.iter){
  celltype.labels <- sce.iter$cellType_broad_hc
  colData(sce.iter)$k2 <- ifelse(grepl("Inhib|Excit", celltype.labels), "neuron",
                                 ifelse(grepl("Oligo|OPC|Astro|EndoMural", celltype.labels), "glial", 
                                        "other"))
  colData(sce.iter)$k3 <- ifelse(grepl("Excit", celltype.labels), "Excit",
                                 ifelse(grepl("Inhib", celltype.labels), "Inhib", 
                                              ifelse(grepl("Oligo|OPC|Astro|EndoMural", celltype.labels), "glial", 
                                                     "other")))
  colData(sce.iter)$k4 <- ifelse(grepl("Excit", celltype.labels), "Excit",
                                 ifelse(grepl("Inhib", celltype.labels), "Inhib", 
                                              ifelse(grepl("Oligo", celltype.labels), "Oligo", 
                                                           ifelse(grepl("OPC|Astro|EndoMural", celltype.labels), "glial_non_oligo", 
                                                                  "other"))))
  return(sce.iter)
})
# name list
names(list.sce.validate.markers) <- c("k2", "k3", "k4")
# save list
save(list.sce.validate.markers, file = file.path(base.path, "list-sce-validate_markers-k-2-3-4_cohort1.rda"))