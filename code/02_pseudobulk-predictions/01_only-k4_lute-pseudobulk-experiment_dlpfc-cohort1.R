#!/usr/bin/env R

#
# Test independent pseudobulk from DLPFC cohort1 snRNAseq data.
#

source("deconvo_method-paper/code/02_pseudobulk-predictions/00_parameters-pseudobulk.R")
sapply(libv, library, character.only = T)
list.sce.markers <- get(load(sce.markers.list.path))
sce <- list.sce.markers$k4

# get experiment results tables
dfp.tall <- get_ypb_experiment_series(sce, sample.id.variable = "Sample", 
                                      celltype.variable = "k4", assay.name = "logcounts",
                                      s.vector = c("Excit" = 10, "Inhib" = 10, "non_oligo_glial" = 3, "Oligo" = 3),
                                      algorithm.name = "nnls", return.dimensions = "tall")



sce
sample.id.variable = "Sample"
celltype.variable = "k4"
assay.name = "logcounts"
s.vector = c("Excit" = 10, "Inhib" = 10, "non_oligo_glial" = 3, "Oligo" = 3)
algorithm.name = "nnls"
return.dimensions = "tall"


# define experiment function
get_ypb_experiment_series <- function(sce, sample.id.variable = "Sample", 
                                      celltype.variable = "k2", assay.name = "logcounts",
                                      s.vector = c("glial" = 3, "neuron" = 10),
                                      algorithm.name = "nnls", return.dimensions = c("wide", "tall"),
                                      dfp.tall.errors = TRUE){
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
    if(dfp.tall.errors){dfp <- dfp_tall_append_errors(dfp)}
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










dfp.wide <- get_ypb_experiment_series(sce, sample.id.variable = "Sample", 
                                      celltype.variable = "k3", assay.name = "logcounts",
                                      s.vector = c("glial" = 3, "Excit" = 10, "Inhib" = 10),
                                      algorithm.name = "nnls", return.dimensions = "wide")

dfp.ct <- dfp_tall_by_celltype(dfp.wide) # dfp.wide, tall by cell type

# make new plots
# plot proportions multipanel
ggplot(dfp.ct, aes(x = true.noscale, y = pred.noscale)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_hline(yintercept = 0.5) + geom_vline(xintercept = 0.5) + theme_bw() +
  xlab("True proportion") + ylab("Predicted proportion") +
  xlim(0, 1) + ylim(0, 1) + facet_wrap(~celltype)
