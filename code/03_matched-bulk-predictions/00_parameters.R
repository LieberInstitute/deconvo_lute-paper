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

# 01 helper functions

pseudobulk_from_sce <- function(sce, group.variable = "donor", 
                            cell.type.variable = "k2", 
                            s = c("glial" = 3, "neuron" = 10), 
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

signature_matrix_from_sce <- function(sce, 
                                      cell.type.variable = "k2", 
                                      summary.method = "mean", 
                                      assay.name = "logcounts"){
  # gets the z signature matrix from an sce object
  lc.condition <- assay.name == "logcounts" & !assay.name %in% names(assays(sce))
  if(lc.condition){sce <- scuttle::logNormCounts(sce)}
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

run_pseudobulk_experiment <- function(list.pb, s, z, method = "nnlsParam"){
  # run an independent psedobulk experiment.
  lapply(list.pb, function(data){
    param.string <- paste0(method, "(y = data$y, z = z, s = s)")
    deconvolution.parameters <- eval(parse(text = param.string))
    p.predicted.pre <- deconvolution.parameters %>% deconvolution()
    p.predicted.final <- p.predicted.pre/sum(p.predicted.pre)
    if(!identical(colnames(p.predicted.final), names(data$p))){
      stop("Error, couldn't match cell type names in data$p, p.predicted.final.")}
    bias <- as.vector(p.predicted.final) - data$p
    rmse <- sqrt(sum(bias^2)/length(bias))
    return(list(p.predicted = p.predicted.final, p.true = data$p, 
                bias = bias, rmse = rmse, 
                parameters = deconvolution.parameters))
  })
}

# 03 helper functions
cell_size_sample <- function(image.table, sample.id, image.sample.id.vector){
  # get k2 cell sizes for a single sample/slide from image outputs
  filter.sample <- image.sample.id.vector==sample.id
  image.table.sample <- image.table[filter.sample,]
  # get k2 labels
  cell.type.vector <- image.table.sample$cell_type
  image.table.sample$k2.type <- ifelse(cell.type.vector %in% c("Inhib", "Excit"), "neuron",
                                       ifelse(cell.type.vector %in% 
                                                c("Oligo", "OPC", "Astro", "Micro"), 
                                              "glial", "other"))
  image.table.sample <- image.table.sample[!image.table.sample$k2.type=="other",]
  # get nuclear area
  by.list <- list(k2.type = image.table.sample$k2.type)
  sample.area <- aggregate(image.table.sample$Nucleus_Area, by = by.list, FUN = "mean")
  sample.area <- c("glial" = sample.area$x[1], "neuron" = sample.area$x[2])
  # get akt3 copies
  by.list <- list(k2.type = image.table.sample$k2.type)
  sample.counts <- aggregate(image.table.sample$AKT3_Copies, by = by.list, FUN = "mean")
  sample.counts <- c("glial" = sample.counts$x[1], "neuron" = sample.counts$x[2])
  return(list(nucleus.area = sample.area, akt3.counts = sample.counts))
}

# 04 helper functions
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

get_k2_area <- function(sizes.table, area.variable, 
                        cell.size.variable = "cell.type", summary.function = "mean"){
  sizes.table$k2 <- ifelse(
    sizes.table[, cell.size.variable] %in% c("Excit", "Inhib"), "neuron",
    ifelse(sizes.table[, cell.size.variable] %in% c("Micro", "Oligo", "Astro"),
           "glial", "other"))
  k2.table <- aggregate(sizes.table[, area.variable], 
                        by = list(k2 = sizes.table$k2), 
                        FUN = summary.function)
  colnames(k2.table) <- c("k2.cell.type", 
                          paste0(area.variable, ".", summary.function))
  return(k2.table)
}

# 05 helper functions
run_deconvolution <- function(y, z, s, method.string){
  param.text <- paste0(method.string, "(y = y, z = z, s = s)")
  new.parameters <- eval(parse(text = param.text))
  predicted.proportions <- deconvolution(new.parameters) %>% suppressMessages()
  predicted.proportions <- predicted.proportions/sum(predicted.proportions)
  return(predicted.proportions)
}

results_series_table <- function(y, z, method.vector, sizes.list){
  # map method and sizes
  table.map <- data.frame(method = rep(method.vector, length(sizes.list)),
                          size = rep(names(sizes.list), 
                                     each = length(method.vector)))
  map.index.vector <- seq(nrow(table.map))
  # get results table
  results.mapped <- lapply(map.index.vector, function(map.index){
    y.index.vector <- y %>% ncol() %>% seq()
    cell.sizes.label <- table.map[map.index,2]
    method.string <- table.map[map.index,1]
    result.iter <- lapply(y.index.vector, function(y.index){
      s <- sizes.list[[cell.sizes.label]]
      y.iter <- y[,y.index,drop=FALSE]
      run_deconvolution(y = y.iter, 
                        z = z.main,
                        s = s, 
                        method.string = method.string)
    })
    result.iter <- do.call(rbind, result.iter) %>% 
      as.data.frame()
    colnames(result.iter)[1:2] <- paste0(colnames(result.iter)[1:2],
                                         ".proportion.predicted")
    result.iter$sample.id <- colnames(y)
    result.iter$cell.sizes.label <- cell.sizes.label
    result.iter$method.string <- method.string
    return(result.iter)
  })
  results.mapped <- do.call(rbind, results.mapped) %>% as.data.frame()
  return(results.mapped)
}

# load cell sizes
cell.sizes.k2.path <- here("deconvo_method-paper", "outputs", "07_cell-size-estimates")

cell.sizes.k2.path <- here(cell.sizes.k2.path, "cell-sizes-k2-table.rda")

# deconvolution algorithms/methods to test
method.vector <- c("nnlsParam", "musicParam", "deconrnaseqParam")

#------------------
# script parameters
#------------------
# 01, independent pseudobulk
# new dlpfc markers
sce.markers.filename <- "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda"
sce.markers.list.path <- here("deconvo_method-paper", "outputs", "09_manuscript", sce.markers.filename)
# mrb dlpfc markers
sce.mrb.name <- "sce-mrb_dlpfc.rda"
# sce.mrb.path <- here("deconvo_method-paper", "outputs", "09_manuscript", sce.mrb.name)
sce.mrb.path <- here("deconvo_method-paper", "outputs", "01_prepare-datasets", sce.mrb.name)

# save results table
independent.pb.results.table.name <- "results-table_independent-pb-mrb.rda"
independent.pb.results.table.path <- here(save.path, independent.pb.results.table.name)

# 02, plot independent pseudobulk results
# save new plots
# scatterplot, proportions
pb.scatterplot.proportions.bysample.colmethod.name <- 
  "ggplot-scatter_neuron-proportions_bysample-colmethod_independent-pb.jpg"
pb.scatterplot.proportions.bysample.colmethod.path <- here(save.path, pb.scatterplot.proportions.bysample.colmethod.name)
# barplot, abs error
pb.barplot.abserror.bysample.xmethod.name <- 
  "ggplot-barplot_neuron-abs-error_bysample-colmethod_independent-pb.jpg"
pb.barplot.abserror.bysample.colmethod.path <- here(save.path, pb.barplot.abserror.bysample.xmethod.name)
# barplot, rmse
pb.barplot.rmse.bysample.xmethod.name <- 
  "ggplot-barplot_rmse_bysample-colmethod_independent-pb.jpg"
pb.barplot.rmse.bysample.colmethod.path <- here(save.path, pb.barplot.rmse.bysample.xmethod.name)

# 03, get image cell proportions
# set the halo data path
halo.output.path <- here("Human_DLPFC_Deconvolution", "processed-data", "03_HALO", "halo_all.Rdata")
# bulk data
rse.k2markers.filepath <- here("deconvo_method-paper", "outputs", "11_bulk-expression-analysis", "rse_k2-marker-expression_ro1-dlpfc.rda")
# data for experiments
assay.name.rse <- "logcounts"
# image cells path
image.cells.name <- "image-cell-counts-table_by-brnum.rda"
image.cells.path <- here(save.path, image.cells.name)
# lexperiment list object
lexperiment.withinsample.name <- "list-experiment-info_within-samples.rda"
lexperiment.withinsample.path <- here(save.path, lexperiment.withinsample.name)
# results table
within.samples.results.table.name <- "results-table_within-sample-deconvolution.rda"
within.samples.results.table.path <- here(save.path, within.samples.results.table.name)
#
plot.barplot.absdiff.byslide.path <- here(save.path, "ggplot-barplot_neuron-prop-abs-diff_by-slide.jpg")
plot.jitterbox.diff.byposition.path <- here(save.path, "ggplot-jitterbox_neuron-prop-diff_by-position.jpg")
plot.jitterbox.absdiff.byposition.path <- here(save.path, "ggplot-jitterbox_neuron-prop-abs-diff_by-position.jpg")
plot.jitterbox.diff.bydonor.path <- here(save.path, "ggplot-jitterbox_neuron-prop-diff_by-donor-brnum.jpg")
plot.jitterbox.absdiff.bydonor.path <- here(save.path, "ggplot-jitterbox_neuron-prop-abs-diff_by-donor-brnum.jpg")
plot.jitterbox.diff.all.path <- here(save.path, "ggplot-jitterbox_neuron-prop-diff_all.jpg")
plot.jitterbox.absdiff.all.path <- here(save.path, "ggplot-jitterbox_neuron-prop-abs-diff_all.jpg")
# load cell sizes list
image.cell.sizes.save.name <- "image_cell-sizes.rda"
image.cell.sizes.save.path <- here("deconvo_method-paper", "outputs", "07_cell-size-estimates", 
                                   image.cell.sizes.save.name)

# 04, get experiment info for within-samples experiments

# 05

# 06, plot within-sample results
# scatterplots
scatterplot.proportions.bysample.colmethod.name <- 
  "ggplot-scatterplot_neuron-proportions_bysample-colmethod.jpg"
scatterplot.proportions.bysample.colmethod.path <- here(save.path, scatterplot.proportions.bysample.colmethod.name)
#
scatterplot.proportions.bylibprep.colmethod.name <- 
  "ggplot-scatterplot_neuron-proportions_bylibprep-colmethod.jpg"
scatterplot.proportions.bylibprep.colmethod.path <- here(save.path, scatterplot.proportions.bylibprep.colmethod.name)
#
scatterplot.proportions.bylibtype.colmethod.name <- 
  "ggplot-scatterplot_neuron-proportions_bylibtype-colmethod.jpg"
scatterplot.proportions.bylibtype.colmethod.path <- here(save.path, scatterplot.proportions.bylibtype.colmethod.name)
#
scatterplot.proportions.byexptgroup.colmethod.name <-
  "ggplot-scatterplot_neuron-proportions_byexptgroup-colmethod.jpg"
scatterplot.proportions.byexptgroup.colmethod.path <- here(save.path, scatterplot.proportions.byexptgroup.colmethod.name)
#
scatterplot.proportions.bymethod.colmethod.name <-
  "ggplot-scatterplot_neuron-proportions_bymethod-colmethod.jpg"
scatterplot.proportions.bymethod.colmethod.path <- here(save.path, scatterplot.proportions.byexptgroup.colmethod.name)

# jitter-boxplots
jitterbox.abserror.bysample.xmethod.name <- 
  "ggplot-jitter-boxplot_neuron-abs-error_bysample-xmethod.jpg"
jitterbox.abserror.bysample.xmethod.path <- here(save.path, jitterbox.abserror.bysample.xmethod.name)
#
jitterbox.abserror.bylibprep.xmethod.name <- 
  "ggplot-jitter-boxplot_neuron-abs-error_bylibprep-xmethod.jpg"
jitterbox.abserror.bylibprep.xmethod.path <- here(save.path, jitterbox.abserror.bylibprep.xmethod.name)
#
jitterbox.abserror.bylibtype.xmethod.name <-
  "ggplot-jitter-boxplot_neuron-abs-error_bylibtype-xmethod.jpg"
jitterbox.abserror.bylibtype.xmethod.path <- here(save.path, jitterbox.abserror.bylibtype.xmethod.name)
#
jitterbox.abserror.byexptgroup.xmethod.name <-
  "ggplot-jitter-boxplot_neuron-abs-error_byexptgroup-xmethod.jpg"
jitterbox.abserror.byexptgroup.xmethod.path <- here(save.path, jitterbox.abserror.byexptgroup.xmethod.name)
#
jitterbox.abserror.bysizetype.xmethod.name <-
  "ggplot-jitter-boxplot_neuron-abs-error_bysizetype-xmethod.jpg"
jitterbox.abserror.bysizetype.xmethod.path <- here(save.path, jitterbox.abserror.bysizetype.xmethod.name)
#
jitterbox.rmse.bysample.xmethod.name <-
  "ggplot-jitter-boxplot_rmse-k2_bysample-xmethod.jpg"
jitterbox.rmse.bysample.xmethod.path <- here(save.path, jitterbox.rmse.bysample.xmethod.name)
#
jitterbox.rmse.bylibprep.xmethod.name <-
  "ggplot-jitter-boxplot_rmse-k2_bylibprep-xmethod.jpg"
jitterbox.rmse.bylibprep.xmethod.path <- here(save.path, jitterbox.rmse.bylibprep.xmethod.name)
#
jitterbox.rmse.bylibtype.xmethod.name <-
  "ggplot-jitter-boxplot_rmse-k2_bylibtype-xmethod.jpg"
jitterbox.rmse.bylibtype.xmethod.path <- here(save.path, jitterbox.rmse.bylibtype.xmethod.name)
#
jitterbox.rmse.byexptgroup.xmethod.name <-
  "ggplot-jitter-boxplot_rmse-k2_byexptgroup-xmethod.jpg"
jitterbox.rmse.byexptgroup.xmethod.path <- here(save.path, jitterbox.rmse.byexptgroup.xmethod.name)

# 04
complete.sample.id.vector.name <- "complete-sample-id-vector_within-samples-trials.rda"
complete.sample.id.vector.path <- here(save.path, complete.sample.id.vector.name)