#!/usr/bin/env R

# Author: Sean Maden
#
# Get cell sizes from preprocessed snRNAseq data.

source("deconvo_method-paper/code/07_cell-size-estimates/00_parameters.R")
sapply(libv, library, character.only = T)
sce <- get(load(sce.path))
# get cell size estimates
sce_summary <- function(sce, expression.summary.type = c("total.counts", "expressed.genes"), 
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
      summary.vector <- colSums(expression.filtered) 
      summary.variable.name <- paste0("total.expression.",assay.name)
    } else{
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
list.total.expression.data <- sce_summary(sce, expression.summary.type = "total.counts")
list.expressed.genes.data <- sce_summary(sce, expression.summary.type = "expressed.genes")
# save
sce.sizes.total.expression.name <- "list-sce-sizes_total-expression.rda"
sce.sizes.total.expression.path <- file.path("deconvo_method-paper", "outputs", 
                                             "07_cell-size-estimates",
                                             sce.sizes.total.expression.name)
sce.sizes.expressed.genes.name <- "list-sce-sizes_expressed-genes.rda"
sce.sizes.expressed.genes.path <- file.path("deconvo_method-paper", "outputs", 
                                            "07_cell-size-estimates",
                                            sce.sizes.expressed.genes.name)
save(list.total.expression.data, file = sce.sizes.total.expression.path)
save(list.expressed.genes.data, file = sce.sizes.expressed.genes.path)