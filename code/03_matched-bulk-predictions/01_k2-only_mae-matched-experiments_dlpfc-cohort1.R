#!/usr/bin/env R

#
#
#

source("deconvo_method-paper/code/03_matched-bulk-predictions/00_parameters-matched-assays.R")
sapply(libv, library, character.only = T)

# load mae 
mae.filepath <- here("deconvo_method-paper", "outputs", "01_prepare-datasets", "mae_final.rda")
mae <- get(load(mae.filepath))
# unpack mae
rse.all <- mae[[1]]
rownames(rse.all) <- rowData(rse.all)$Symbol
# snrnaseq reference -- using same reference across experiments
sce.iter <- mae[[2]]
sce.iter <- logNormCounts(sce.iter)
# experiment variables
assay.name <- "logcounts"
deconvolution.algorithm <- "nnls"

# get cell sizes
# get cell size factor series
# load cell size scale factors
df.csf <- get_csf_reference()
df.csf <- df.csf[df.csf$tissue=="brain",]
df.csf.area <- df.csf[grepl("mRNA", df.csf$scale.factor.type),]
df.csf.mrna <- df.csf[df.csf$scale.factor.type=="cell area",]
s.set.osm.area <- c("glial" = df.csf.area[df.csf.area$cell_type=="glial",]$scale.factor.value,
                    "neuron" = df.csf.area[df.csf.area$cell_type=="neuron",]$scale.factor.value)
s.set.osm.mrna <- c("glial" = df.csf.mrna[df.csf.mrna$cell_type=="glial",]$scale.factor.value,
                    "neuron" = df.csf.mrna[df.csf.mrna$cell_type=="neuron",]$scale.factor.value)
# coerce to list for experiment
list.s.pred <- list(s.set1 = c("glial" = 3, "neuron" = 10),
                    s.set2 = c("glial" = 10, "neuron" = 3),
                    s.set3 = c("glial" = 1, "neuron" = 1),
                    s.set4 = s.set.osm.area, s.set5 = s.set.osm.mrna)

#---------------------------------------------------
# k2 experiment -- same reference across experiments
#---------------------------------------------------
# run k2 experiment
celltype.variable <- "k2"
sample.id.vector <- unique(colData(mae)[,1][complete.cases(mae)])
sample.id <- sample.id.vector[1]
df.s.k2.shared <- do.call(rbind, lapply(seq(length(list.s.pred)), function(s.index){
  # format cell sizes
  s.set.name <- names(list.s.pred)[s.index]
  s.vector.pred <- list.s.pred[[s.index]]
  s.vector.pred <- order_svector(s.vector.pred)
  # iterate on samples, returning predictions matrix
  do.call(rbind, lapply(sample.id.vector, function(sample.id){
    message(sample.id)
    # filter bulk
    filter.y.sample <- colData(rse.all)[,"batch.id2"]==sample.id
    filter.y.marker <- rownames(rse.all) %in% rownames(sce.iter)
    rse.iter <- rse.all[filter.y.marker, filter.y.sample]
    rse.iter <- logNormCounts(rse.iter)
    y.iter <- assays(rse.iter)[[assay.name]]
    # get predictions
    prop.pred.iter <- lute(sce = sce.iter, y = y.iter, assay.name = assay.name, 
                           celltype.variable = celltype.variable, s = s.vector.pred, 
                           typemarker.algorithm = NULL, return.info = FALSE,
                           deconvolution.algorithm = deconvolution.algorithm)$deconvolution.results %>%
      as.data.frame()
    prop.pred.iter$s.set.label <- s.set.name
    prop.pred.iter$s.set.values <- paste0(
      names(s.vector.pred),"=",s.vector.pred,collapse=",")
    prop.pred.iter$sample.id <- sample.id
    prop.pred.iter$k.type <- celltype.variable
    prop.pred.iter
  }))
}))
df.s.k2.shared$experiment.type <- "shared.snrnaseq.reference"

#--------------------------------------------------------------------------
# k2 experiment -- using within-sample matched reference across experiments
#--------------------------------------------------------------------------
# run k2 experiment
celltype.variable <- "k2"
sample.id.vector <- unique(colData(mae)[,1][complete.cases(mae)])
sample.id <- sample.id.vector[1]
df.s.k2.within <- do.call(rbind, lapply(seq(length(list.s.pred)), function(s.index){
  # format cell sizes
  s.set.name <- names(list.s.pred)[s.index]
  s.vector.pred <- list.s.pred[[s.index]]
  s.vector.pred <- order_svector(s.vector.pred)
  # iterate on samples, returning predictions matrix
  do.call(rbind, lapply(sample.id.vector, function(sample.id){
    message(sample.id)
    # filter bulk
    filter.y.sample <- colData(rse.all)[,"batch.id2"]==sample.id
    filter.y.marker <- rownames(rse.all) %in% rownames(sce.iter)
    rse.iter <- rse.all[filter.y.marker, filter.y.sample]
    rse.iter <- logNormCounts(rse.iter)
    y.iter <- assays(rse.iter)[[assay.name]]
    # filter snrnaseq
    sce.iter.sample <- sce.iter[,sce.iter$Sample==sample.id]
    # get predictions
    prop.pred.iter <- lute(sce = sce.iter.sample, y = y.iter, assay.name = assay.name, 
                           celltype.variable = celltype.variable, s = s.vector.pred, 
                           typemarker.algorithm = NULL, return.info = FALSE,
                           deconvolution.algorithm = deconvolution.algorithm)$deconvolution.results %>%
      as.data.frame()
    prop.pred.iter$s.set.label <- s.set.name
    prop.pred.iter$s.set.values <- paste0(
      names(s.vector.pred),"=",s.vector.pred,collapse=",")
    prop.pred.iter$sample.id <- sample.id
    prop.pred.iter$k.type <- celltype.variable
    prop.pred.iter
  }))
}))
df.s.k2.within$experiment.type <- "within.sample.snrnaseq.reference"

#---------------
# format results
#---------------
df.k2 <- rbind(df.s.k2.shared, df.s.k2.within)




