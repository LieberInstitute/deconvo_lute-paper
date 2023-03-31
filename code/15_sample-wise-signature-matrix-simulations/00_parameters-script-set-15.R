#!/usr/bin/env R

# Author: Sean Maden
#
# Main parameters, or dependency objects, for sample-wise signatures trials.


# dependencies
libv <- c("here", "lute", "dplyr", "ggplot2", "gridExtra", "SingleCellExperiment", 
          "SummarizedExperiment", "scran", "DeconvoBuddies", "UpSetR")
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
    marker_genes(sce[,filter], celltype.variable, assay.name, genes.per.type)
  })
  names(batch.markers.list) <- unique.batch.vector
  message("Finished all marker lists.")
  return(batch.markers.list)
}

marker_genes <- function(sce, celltype.variable, assay.name, genes.per.type = 20){
  message("Running marker tests...")
  mr <- DeconvoBuddies::get_mean_ratio2(sce, assay_name = assay.name, cellType_col = celltype.variable)
  mr <- mr %>% as.data.frame()
  message("Getting top ",genes.per.type," markers per cell type...")
  unique.cell.types <- mr[,"cellType.target"] %>% as.character() %>% unique()
  top.markers.list <- lapply(unique.cell.types, function(unique.type.id){
    mr %>% filter(cellType.target == unique.type.id) %>% arrange(rank_ratio) %>% 
      top_n(n = genes.per.type)
  })
  top.markers.table <- do.call(rbind, top.markers.list)
  return(top.markers.table)
}

#------------------
# script parameters
#------------------
# 01, independent pseudobulk
# get full sce data
sce.name <- "sce_DLPFC.Rdata"
sce.path <- file.path("DLPFC_snRNAseq/processed-data/sce", sce.name)

# 04 mrb setup
# mrb dlpfc markers
sce.mrb.name <- "sce-mrb_dlpfc.rda"
sce.mrb.path <- here("deconvo_method-paper", "outputs", "09_manuscript", sce.mrb.name)
