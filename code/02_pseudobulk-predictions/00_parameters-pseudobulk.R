#!/usr/bin/env R

# Author: Sean Maden
#
# Main parameters, or dependency objects, for deconvolution framework trials.

# dependencies
libv <- c("here", "lute", "dplyr", "ggplot2", "gridExtra", "SingleCellExperiment", 
          "SummarizedExperiment", "scran")
sapply(libv, library, character.only = TRUE)

# save path
save.path <- here("deconvo_method-paper", "outputs", "02_pseudobulk-predictions")

# mrb sce path
sce.mrb.path <- here("deconvo_method-paper", "outputs", "01_prepare-datasets", "sce-mrb_dlpfc.rda")

# dlpfc markers path
sce.markers.list.path <- here("deconvo_method-paper", "outputs", "02_pseudobulk-predictions", "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda")

# define experiment function
get_ypb_experiment_series <- function(sce, sample.id.variable = "Sample", 
                                      celltype.variable = "k2", assay.name = "logcounts",
                                      s.vector = c("glial" = 3, "neuron" = 10),
                                      algorithm.name = "nnls", return.dimensions = c("wide", "tall")){
  # get pseudobulk experiment series, testing cellsize adjustment
  # use with dfp_tall_by_celltype()
  # get experiment results
  if(!sample.id.variable %in% colnames(colData(sce))){stop("Error: couldn't find sample.id.variable in sce coldata.")}
  celltype.variable.format <- sce[[celltype.variable]] %>% as.factor()
  unique.celltypes <- unique(celltype.variable.format)
  num.celltypes <- length(unique(celltype.variable.format))
  s.vector.null <- rep(1, num.celltypes); names(s.vector.null) <- unique.celltypes
  s.vector.null <- s.vector.null[order(match(names(s.vector.null), names(s.vector)))]
  # get experiment proportions
  dfp.noscale <- get_ypb_experiment_results(sce, 
                                               sample.id.variable = sample.id.variable, 
                                               celltype.variable = celltype.variable, 
                                               assay.name = assay.name,
                                               s.vector.ypb = s.vector,
                                               s.vector.pred = s.vector.null)
  dfp.withscale <- get_ypb_experiment_results(sce, 
                                                 sample.id.variable = sample.id.variable, 
                                                 celltype.variable = celltype.variable, 
                                                 assay.name = assay.name,
                                                 s.vector.ypb = s.vector,
                                                 s.vector.pred = s.vector)
  if(return.dimensions == "tall"){
    # get plot data -- tall
    dfp.noscale$type <- 'noscale'
    dfp.withscale$type <- 'withscale'
    dfp.noscale$sample.id <- rownames(dfp.noscale)
    dfp.withscale$sample.id <- rownames(dfp.withscale)
    dfp <- rbind(dfp.noscale, dfp.withscale)
  } else{
    # get plot data -- wide
    colnames(dfp.noscale) <- paste0(colnames(dfp.noscale), ".noscale")
    colnames(dfp.withscale) <- paste0(colnames(dfp.withscale), ".withscale")
    identical(rownames(dfp.noscale), rownames(dfp.withscale))
    dfp <- cbind(dfp.noscale, dfp.withscale)
  }
  return(dfp)
}

get_ypb_experiment_results <- function(sce, sample.id.variable = "Sample", 
                                       celltype.variable = "k2", assay.name = "logcounts",
                                       s.vector.ypb = c("glial" = 3, "neuron" = 10),
                                       s.vector.pred = c("glial" = 1, "neuron" = 1),
                                       algorithm.name = "nnls"){
  # get results for a single iteration of an experiment
  # use with get_ypb_experiment_series()
  if(assay.name == "logcounts" & !"logcounts" %in% names(assays(sce))){sce <- scuttle::logNormCounts(sce)}
  unique.sample.id.vector <- unique(sce[[sample.id.variable]])
  dfp <- do.call(rbind, lapply(unique.sample.id.vector, function(sample.id){
    sce.iter <- sce[,sce[[sample.id.variable]]==sample.id]
    ypb.iter <- ypb_from_sce(sce = sce.iter, assay.name = assay.name, 
                             celltype.variable = celltype.variable,
                             sample.id.variable = sample.id.variable, 
                             S = s.vector.ypb) %>% as.matrix()
    prop.true.iter <- table(sce.iter[[celltype.variable]]) %>% prop.table() %>% as.matrix() %>% t()
    prop.pred.iter <- lute(sce = sce.iter, y = ypb.iter, assay.name = assay.name, 
                           celltype.variable = celltype.variable, s = s.vector.pred, 
                           typemarker.algorithm = NULL, return.info = FALSE,
                           deconvolution.algorithm = algorithm.name)$deconvolution.results
    colnames(prop.pred.iter) <- paste0(colnames(prop.pred.iter), ".pred")
    colnames(prop.true.iter) <- paste0(colnames(prop.true.iter), ".true")
    dfp.iter <- cbind(prop.true.iter, prop.pred.iter) %>% as.data.frame()
    dfp.iter
  }))
  rownames(dfp) <- unique.sample.id.vector
  return(dfp)
}

dfp_tall_by_celltype <- function(dfp.wide){
  # input: dfp.wide from get_ypb_experiment_series()
  cn.wide <- colnames(dfp.wide)
  ct.wide <- unique(gsub("\\..*", "", cn.wide))
  dfp.tall.by.celltype <- do.call(rbind, lapply(ct.wide, function(ct.iter){
    dfp.wide.iter <- dfp.wide[,grepl(paste0(ct.iter,"\\..*"), cn.wide)]
    colnames(dfp.wide.iter) <- gsub(paste0(ct.iter, "\\."), "", colnames(dfp.wide.iter))
    dfp.wide.iter$celltype <- ct.iter
    dfp.wide.iter$sample.id <- rownames(dfp.wide.iter)
    dfp.wide.iter
  }))
  return(dfp.tall.by.celltype)
}
