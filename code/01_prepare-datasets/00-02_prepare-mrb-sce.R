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
celltypevar <- "cellType"
# define marker categories
sce[["k2"]] <- ifelse(grepl("^Excit.*|^Inhib.*", sce[[celltypevar]]), "neuron", "other")
sce[["k3"]] <- ifelse(grepl("^Excit.*", sce[[celltypevar]]), "Excit", 
                      ifelse(grepl("^Inhib.*", sce[[celltypevar]]), "Inhib", "other"))
sce[["k4"]] <- ifelse(grepl("^Excit.*", sce[[celltypevar]]), "Excit", 
                      ifelse(grepl("^Inhib.*", sce[[celltypevar]]), "Inhib", 
                             ifelse(grepl("^Oligo$", sce[[celltypevar]]), "Oligo", "other")))

#--------------
# save sce data
#--------------
save(sce, file = save.path)