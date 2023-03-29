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

# helper functions

pseudobulk_from_sce <- function(sce, group.variable = "donor", 
                            cell.type.variable = "k2", 
                            s = c("neuron" = 3, "glial" = 10), 
                            summary.method = "mean", 
                            assay.name = "logcounts"){
  # makes a new pseudobulk sample from an sce object
  cd <- colData(sce)
  unique.groups <- unique(cd[,group.variable])
  unique.cell.types <- unique(cd[,cell.type.variable])
  unique.cell.types <- unique.cell.types[order(unique.cell.types)]
  s <- s[order(names(s))]
  if(!identical(names(s), unique.cell.types)){
    stop("Error, couldn't match unique cell types to s names.")}
  lpb <- lapply(unique.groups, function(group.index){
    # group.index <- unique.groups[1]
    filter <- cd[,group.variable] == group.index
    scef <- sce[,filter]; cdf <- colData(scef)
    # get zs transformed signature matrix
    expression.matrix <- assays(scef)[[assay.name]]
    z <- do.call(cbind, lapply(unique.cell.types, function(cell.type.index){
      index.filter <- cdf[,cell.type.variable]==cell.type.index
      rowMeans(expression.matrix[,index.filter])
    }))
    colnames(z) <- unique.cell.types
    zs <- sweep(z, 2, STATS = s, FUN = "*")
    # get p cell type proportions
    p.true.counts <- table(cdf[,cell.type.variable])
    p.true.counts <- p.true.counts[order(names(p.true.counts))]
    if(!identical(names(p.true.counts), colnames(zs))){
      stop("Error, couldn't match cell type labels in p.true and zs")}
    p.true.proportions <- p <- prop.table(p.true.counts)
    # final pseudobulk
    y.pseudobulk <- t(t(p) %*% t(zs))
    return(list(y = y.pseudobulk, p = p.true.proportions, 
                p.counts = p.true.counts))
  })
  names(lpb) <- unique.groups
  return(lpb)
}

signature_matrix_from_sce <- function(sce, cell.type.variable = "k2", 
                                      summary.method = "mean", 
                                      assay.name = "logcounts"){
  # gets the z signature matrix from an sce object
  expression.matrix <- assays(sce)[[assay.name]] %>% as.matrix()
  cd <- colData(sce)
  unique.cell.types <- cd[,cell.type.variable] %>% unique()
  unique.cell.types <- unique.cell.types[order(unique.cell.types)]
  z <- do.call(cbind, lapply(unique.cell.types, function(cell.type.index){
    filter.index <- cd[,cell.type.variable]==cell.type.index
    if(summary.method == "mean"){
      rowMeans(expression.matrix[,filter.index])
    } else{stop('Error, unrecognized summary.method.')}
  }))
  colnames(z) <- unique.cell.types
  return(z)
}

run_experiment <- function(list.pb, method = "nnlsParam"){
  # run an independent psedobulk experiment.
  lapply(list.pb, function(data){
    param.string <- paste0(method, "(y = data$y, z = z, s = s)")
    deconvolution.parameters <- eval(parse(text = param.string))
    p.predicted.pre <- deconvolution.parameters %>% deconvolution()
    p.predicted.final <- p.predicted.pre/sum(p.predicted.pre)
    if(!identical(names(p.predicted.final), names(data$p))){
      stop("Error, couldn't match cell type names in data$p, p.predicted.final.")}
    bias <- p.predicted.final - data$p
    rmse <- sqrt(sum(bias^2)/length(bias))
    return(list(p.predicted = p.predicted.final, p.true = data$p, 
                bias = bias, rmse = rmse, parameters = deconvolution.parameters))
  })
}

standard_sample_id <- function(table, location.label, brnum.label = "BrNum"){
  # get standard formatted sample identifiers from sn, bulk, and image tables.
  brnum.vector <- table[,brnum.label]
  location.vector <- toupper(table[,location.label])
  location.vector <- ifelse(grepl("ANT", location.vector), "ANTERIOR",
                            ifelse(grepl("MID", location.vector), "MIDDLE",
                                   ifelse(grepl("POST", location.vector), "POSTERIOR", 
                                          location.vector)))
  return(paste0(brnum.vector, "_", location.vector))
}

# load cell sizes
cell.sizes.k2.path <- here("deconvo_method-paper", "outputs", "07_cell-size-estimates")
cell.sizes.k2.path <- here(cell.sizes.k2.path, "cell-sizes-k2-table.rda")

#------------------
# script parameters
#------------------
# 01, independent pseudobulk
# new dlpfc markers
sce.markers.filename <- "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda"
sce.markers.list.path <- here("deconvo_method-paper", "outputs", "09_manuscript", sce.markers.filename)
# mrb dlpfc markers
sce.mrb.name <- "sce-mrb_dlpfc.rda"
sce.mrb.path <- here("deconvo_method-paper", "outputs", "09_manuscript", sce.mrb.name)

# 02, within-samples tests
# set the halo data path
halo.output.path <- here("Human_DLPFC_Deconvolution", "processed-data", "03_HALO", "halo_all.Rdata")
# bulk data
rse.k2markers.filepath <- here("deconvo_method-paper", "outputs", "11_bulk-expression-analysis", "rse_k2-marker-expression_ro1-dlpfc.rda")
# data for experiments
assay.name.rse <- "counts"


# 03, across-samples tests

# 04