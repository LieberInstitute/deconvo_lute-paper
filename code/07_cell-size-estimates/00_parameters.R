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

sce_summary <- function(sce, expression.summary.type = c("total.expression", "expressed.genes"), 
                        assay.name = "logcounts", cell.types.variable = "cellType_broad_hc",
                        size.summary.type = "median", batch.id = "SAMPLE_ID",
                        minimum.expression = 2){
  expression.delayed.matrix <- assays(sce)[[assay.name]]
  unique.types <- sce[[cell.types.variable]] %>% unique()
  cd <- sce %>% colData()
  # get expression summary per expression.summary.type
  expression.summary <- do.call(rbind, lapply(unique.types, function(cell.type){
    message("Working on cell type: ", cell.type, "...")
    filter <- cd[,cell.types.variable] == cell.type
    expression.filtered <- expression.delayed.matrix[,filter] %>% as.matrix()
    # parse summary type option
    if(expression.summary.type == "total.expression"){
      message("getting total expression...")
      summary.vector <- colSums(expression.filtered) 
      summary.variable.name <- paste0("total.expression.",assay.name)
    } else{
      message("getting expressed genes...")
      summary.vector <- apply(expression.filtered, 2, function(cell.data){
        length(cell.data[cell.data > minimum.expression])
      })
      summary.variable.name <- paste0("expressed.genes.min.",
                                      minimum.expression,".",assay.name)
    }
    # make summary table
    summary.table <- summary.vector %>% matrix(ncol = 1) %>% as.data.frame()
    summary.table$cell.barcode <- names(summary.vector)
    colnames(summary.table) <- summary.variable.name
    summary.table$cell.type <- cell.type
    summary.table$batch.id <- cd[filter, batch.id]
    summary.table
  }))
  message("finished with expression summary.")
  # get cell size estimates per size.summary.type
  by.list <- list(cell.type = expression.summary$cell.type)
  cell.size.summary <- aggregate(expression.summary[,1], by = by.list, 
                                 FUN = size.summary.type)
  # return all
  return(list(expression.summary.type = expression.summary.type,
              expression.summary = expression.summary,
              cell.size.summary = cell.size.summary))
}

#---------------------
# parameters by script
#---------------------
# 01 image cell sizes
area.variables <- c("Nucleus_Area", "AKT3_Copies") # table columns to summarize
image.cell.sizes.save.name <- "image_cell-sizes.rda"
image.cell.sizes.save.path <- here(save.path, image.cell.sizes.save.name)

# 02 snrnaseq cell sizes
sce.file.name <- "sce_DLPFC.Rdata"
sce.path <- here("DLPFC_snRNAseq","processed-data", "sce", sce.file.name)
# save paths
sce.sizes.total.expression.name <- "list-sce-sizes_total-expression.rda"
sce.sizes.total.expression.path <- file.path("deconvo_method-paper", "outputs", 
                                             "07_cell-size-estimates",
                                             sce.sizes.total.expression.name)
sce.sizes.expressed.genes.name <- "list-sce-sizes_expressed-genes.rda"
sce.sizes.expressed.genes.path <- file.path("deconvo_method-paper", "outputs", 
                                            "07_cell-size-estimates",
                                            sce.sizes.expressed.genes.name)

# 03
