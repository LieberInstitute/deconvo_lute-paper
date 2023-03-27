#!/usr/bin/env R

# Author: Sean Maden
#
# Parameters for calculation and cell size estimates from snRNAseq, image analysis.
#

libv <- c("here", "dplyr", "ggforce", "ggplot2", "gridExtra", "ggpubr", 
          "data.table", "ggcorrplot", "SingleCellExperiment", 
          "SummarizedExperiment")
sapply(libv, library, character.only = TRUE)

# set the save directory
save.path <- here("deconvo_method-paper", "outputs", "07_cell-size-estimates")

# set the halo data path
halo.output.file.name <- "halo_all.Rdata"
halo.output.path <- here("Human_DLPFC_Deconvolution", "processed-data", "03_HALO", 
                         halo.output.file.name)

# cell labels
labels <- c("Endo" = "CLDN5", "Astro" = "GFAP", "Inhib" = "GAD1", 
            "Excit" = "SLC17A7", "Micro" = "TMEM119", "Oligo" = "OLIG2")

#-----------------
# helper functions
#-----------------
cell_sizes <- function(data, area.variable = "Nucleus_Area", 
                       by.variable = "cell_type", fun = "median"){
  aggregate(data[,area.variable], by = list(variable = data[,by.variable]), FUN = fun)
}

#---------------------
# parameters by script
#---------------------
# 01 image cell sizes
image.cell.sizes.save.name <- "image_cell-sizes.rda"
image.cell.sizes.save.path <- here(save.path, image.cell.sizes.save.name)

# 02 snrnaseq cell sizes
sce.file.name <- "sce_DLPFC.Rdata"
sce.path <- here("DLPFC_snRNAseq","processed-data", "sce", sce.file.name)

sn.cell.sizes.save.name <- "snrnaseq_cell-sizes.rda"
sn.cell.sizes.save.path <- here(save.path, sn.cell.sizes.save.name)

