#!/usr/bin/env R

# Author: Sean Maden
#
# Main parameters, or dependency objects, for deconvolution framework trials.

# dependencies
libv <- c("here", "lute", "dplyr", "ggplot2", "gridExtra", "SingleCellExperiment", 
          "SummarizedExperiment", "scran")
sapply(libv, library, character.only = TRUE)

# save path
save.path <- here("deconvo_method-paper", "outputs", "14_deconvolution-framework-trials")

#-----------------
# helper functions
#-----------------

append_k_columns <- function(df.input, df.ct = NULL, celltype.variable = "cell_type"){
  df.ct <- data.frame(cell_type = c("Inhib", "Other", "Astro", "Endo", "Excit", "Oligo", "OPC", "Micro", "all"),
                      k2 = c("neuron", "NA", "glial", "glial", "neuron", "glial", "glial", "glial", "NA"),
                      k3 = c("Inhib", "NA", "glial", "glial", "Excit", "glial", "glial", "glial", "NA"),
                      k4 = c("Inhib", "NA", "non_oligo_glial", "non_oligo_glial", "Excit", "Oligo", 
                             "non_oligo_glial", "non_oligo_glial", "NA"))
  celltype.vector <- df.input[,celltype.variable]
  unique.cell.types <- unique(celltype.vector)
  df.ct.filter <- df.ct[df.ct$cell_type %in% unique.cell.types,]
  df.ct.new <- do.call(rbind, lapply(celltype.vector, function(cell.type.id){
    df.ct.filter[df.ct.filter[,"cell_type"]==cell.type.id,,drop = F]
  }))
  cbind(df.input, df.ct.new)
}

get_ymatch_experiment_results <- function(mae, sample.id.variable = "Sample", 
                                       celltype.variable = "k2", assay.name = "logcounts",
                                       s.vector.pred = c("glial" = 3, "neuron" = 10),
                                       deconvolution.algorithm = "nnls", system.sleep.sec = 2){
  sce <- mae[[2]]
  rse.all <- mae[[1]]
  
  s.vector.ypb <- order_svector(s.vector.ypb)
  s.vector.pred <- order_svector(s.vector.pred)
  # get results for a single iteration of an experiment
  # use with get_ypb_experiment_series()
  if(assay.name == "logcounts" & !"logcounts" %in% names(assays(sce))){sce <- scuttle::logNormCounts(sce)}
  unique.sample.id.vector <- unique(sce[[sample.id.variable]])
  dfp <- do.call(rbind, lapply(unique.sample.id.vector, function(sample.id){
    sce.iter <- sce[,sce[[sample.id.variable]]==sample.id]
    y.iter <- mae[[1]] %>% as.matrix()
    prop.true.iter <- table(sce.iter[[celltype.variable]]) %>% prop.table() %>% as.matrix() %>% t()
    prop.pred.iter <- lute(sce = sce.iter, y = ypb.iter, assay.name = assay.name, 
                           celltype.variable = celltype.variable, s = s.vector.pred, 
                           typemarker.algorithm = NULL, return.info = FALSE,
                           deconvolution.algorithm = deconvolution.algorithm)$deconvolution.results
    colnames(prop.pred.iter) <- paste0(colnames(prop.pred.iter), ".pred")
    colnames(prop.true.iter) <- paste0(colnames(prop.true.iter), ".true")
    dfp.iter <- cbind(prop.true.iter, prop.pred.iter) %>% as.data.frame()
    Sys.sleep(system.sleep.sec)
    dfp.iter
  }))
  rownames(dfp) <- unique.sample.id.vector
  return(dfp)
}

order_svector <- function(s.vector = c("neuron" = 10, "glial" = 3)){
  s.vector[names(s.vector)[order(names(s.vector))]]
}
