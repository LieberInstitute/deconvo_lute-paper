#!/usr/bin/env R

# Author: Sean Maden
#
# Main parameters, or dependency objects, for sample-wise signatures trials.


# dependencies
libv <- c("here", "lute", "dplyr", "ggplot2", "gridExtra", "SingleCellExperiment", 
          "SummarizedExperiment", "scran", "DeconvoBuddies", "UpSetR", "DelayedArray")
sapply(libv, library, character.only = TRUE)

# save path
save.path <- here("deconvo_method-paper", "outputs", 
                  "15_sample-wise-signature-matrix-simulations")


#-----------------
# helper functions
#-----------------
# get marker gene lists
markers_by_batch <- function(sce, batch.variable, celltype.variable, assay.name, genes.per.type){
  unique.batch.vector <- sce[[batch.variable]] %>% unique() %>% as.character()
  # get gene markers for each specified batch
  batch.markers.list <- lapply(unique.batch.vector, function(batch.id){
    message("Getting markers for batch id: ", batch.id, "...")
    filter <- sce[[batch.variable]] == batch.id
    lute(sce = sce[,filter], 
         celltype.variable = celltype.variable, 
         assay.name = assay.name, 
         markers.per.type = genes.per.type,
         deconvolution.algorithm = NULL,
         return.info = TRUE)
  })
  names(batch.markers.list) <- unique.batch.vector
  message("Finished all marker lists.")
  return(batch.markers.list)
}

get.overlapping.markers <- function(markers.by.batch, min.overlap.rate = 0.8){
  # parse cell types
  batch.id.vector <- names(markers.by.batch)
  marker.list <- lapply(batch.id.vector, function(batch.id){
    result.batch <- markers.by.batch[[batch.id]][[1]]
    table.batch <- result.batch$result.info
    table.batch <- table.batch %>% filter(gene %in% result.batch$markers)
    table.batch$batch.id <- batch.id
    return(table.batch)
  })
  marker.table <- do.call(rbind, marker.list) %>% as.data.frame()
  unique.cell.types <- marker.table$cellType.target %>% unique()
  # get overlaps by cell type
  list.markers.final <- lapply(unique.cell.types, function(type){
    marker.table.type <- marker.table %>% filter(cellType.target == type)
    overlap.frequency <- marker.table.type$gene %>% table() %>% as.data.frame()
    overlap.frequency$rate.overlap <- 100*overlap.frequency[,2]/length(markers.by.batch)
    if(is(min.overlap.rate, "NULL")){
      overlap.frequency
    } else{
      filter.overlaps <- overlap.frequency$percent.overlap >= min.overlap.rate
      overlap.frequency[filter.overlaps, 1] %>% as.character()
    }
  })
  names(list.markers.final) <- unique.cell.types
  return(list.markers.final)
}

signature_matrix_from_sce <- function(sce, 
                                      cell.type.variable = "k2", 
                                      summary.method = "mean", 
                                      assay.name = "counts"){
  # gets the z signature matrix from an sce object
  expression.matrix <- assays(sce)[[assay.name]] %>% as.matrix()
  cd <- colData(sce)
  unique.cell.types <- cd[,cell.type.variable] %>% unique()
  unique.cell.types <- unique.cell.types[order(unique.cell.types)]
  z <- do.call(cbind, lapply(unique.cell.types, function(cell.type.index){
    filter.index <- cd[,cell.type.variable]==cell.type.index
    if(summary.method == "mean"){
      DelayedArray::rowMeans(expression.matrix[,filter.index])
    } else{stop('Error, unrecognized summary.method.')}
  }))
  colnames(z) <- unique.cell.types
  return(z)
}

#------------------
# script parameters
#------------------
# 01, independent pseudobulk
# get full sce data
sce.name <- "sce_DLPFC.Rdata"
sce.path <- file.path("DLPFC_snRNAseq/processed-data/sce", sce.name)
sce.prepared.path <- file.path("deconvo_method-paper", "outputs", 
                               "15_sample-wise-signature-matrix-simulations",
                               "sce-prepared_dlpfc-ro1-train.rda")

# 02, markers by slide, overlaps
list.markers.final.name <- "list-markers-by-slide-overlaps_dlpfc-ro1-train.rda"
list.markers.final.path <- here(save.path, list.markers.final.name)

# 04 mrb setup
# mrb dlpfc markers
sce.mrb.name <- "sce-mrb_dlpfc.rda"
sce.mrb.path <- here("deconvo_method-paper", "outputs", "09_manuscript", sce.mrb.name)
