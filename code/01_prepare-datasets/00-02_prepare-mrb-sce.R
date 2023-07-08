#!/usr/bin/env R

#
# Prepares Tran et al 2021 DLPFC snRNAseq data for analyses. Steps:
# * Assigns k cell type labels
#

libv <- c("lute", "scuttle", "dplyr", "limma", "ggplot2", "ggforce", "gridExtra",
          "glmGamPoi", "sva", "DeconvoBuddies", "SingleCellExperiment", "limma",
          "SummarizedExperiment")
sapply(libv, library, character.only = TRUE)

#----------
# load data
#----------
# get save dpath
code.folder <- "01_prepare-datasets"
project.folder <- "deconvo_method-paper"
save.filename <- "sce-mrb_dlpfc.rda"
save.path <- here(project.folder, "outputs", code.folder, save.filename)

# load downloaded sce data (see ReadMe)
sce.downloaded.filepath <- here(project.folder, "data", "sce", "SCE_DLPFC-n3_tran-etal.rda")
sce <- get(load(sce.downloaded.filepath))

#-----------------------------------
# assign marker labels at variable k
#-----------------------------------
# filter non-neuron, non-glial
filter.types.neuron.glial <- grepl("Excit|Inhib|Oligo|OPC|Micro|Astro", sce[["cellType"]])
sce <- sce[,filter.types.neuron.glial]

# define marker categories from cell types vector
mrb.cell.type.vector <- sce[["cellType"]]
mrb.is.neuron <- grepl("Excit|Inhib", mrb.cell.type.vector)
mrb.is.glial <- grepl("Oligo|OPC|Micro|Astro", mrb.cell.type.vector)
mrb.is.oligo <- grepl("Oligo", mrb.cell.type.vector)
mrb.is.non.oligo.glial <- grepl("OPC|Micro|Astro", mrb.cell.type.vector)

mrb.is.excit <- grepl("Excit", mrb.cell.type.vector)
mrb.is.inhib <- grepl("Inhib", mrb.cell.type.vector)

# define k labels
colData(sce)[,"k2"] <- ifelse(mrb.is.neuron, "neuron", 
                                  ifelse(mrb.is.glial, "glial", "other"))
colData(sce)[,"k3"] <- ifelse(mrb.is.excit, "excit", 
                              ifelse(mrb.is.inhib, "inhib", 
                                     ifelse(mrb.is.glial, "glial", "other")))
colData(sce)[,"k4"] <- ifelse(mrb.is.excit, "excit", 
                              ifelse(mrb.is.inhib, "inhib", 
                                     ifelse(mrb.is.oligo, "oligo",
                                          ifelse(mrb.is.non.oligo.glial, "non_oligo_glial", "other"))))

#--------------
# save sce data
#--------------
save(sce, file = save.path)