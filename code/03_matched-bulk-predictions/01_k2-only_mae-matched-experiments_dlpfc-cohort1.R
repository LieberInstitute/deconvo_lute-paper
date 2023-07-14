#!/usr/bin/env R

#
#
#

source("deconvo_method-paper/code/03_matched-bulk-predictions/00_parameters-matched-assays.R")
sapply(libv, library, character.only = T)

# load mae 
# mae.filepath <- here("deconvo_method-paper", "outputs", "01_prepare-datasets", "mae_final.rda")
mae.final.filepath <- here("deconvo_method-paper", "outputs", "01_prepare-datasets", "mae_additional-data_final.rda")
mae <- get(load(mae.final.filepath))
# unpack mae
rse.all <- mae[["bulk.rnaseq"]]
rownames(rse.all) <- rowData(rse.all)$Symbol
# snrnaseq reference -- using same reference across experiments
sce.iter <- mae[["sn1.rnaseq"]]
sce.iter <- logNormCounts(sce.iter)
# experiment variables
assay.name <- "logcounts"
deconvolution.algorithm <- "nnls"
# get true proportions from rnascope data
rnascope <- mae[["rnascope.image"]]
# append k2 label
cd <- colData(rnascope)
cd$k2 <- "NA"
cd$k2 <- ifelse(grepl("Excit|Inhib", cd$cell_type), "neuron",
                ifelse(grepl("Endo|Oligo|Micro", cd$cell_type), "glial", "NA"))
colData(rnascope) <- cd
# filter na
rnascope <- rnascope[,!rnascope$k2=="NA"]
sample.id.vector <- unique(rnascope$Sample)
cell.type.variable <- "k2"
df.rn <- do.call(rbind, lapply(sample.id.vector, function(sample.id){
  rnf <- rnascope[,rnascope$Sample==sample.id]
  # proportions
  df.prop <- table(rnf[[cell.type.variable]], rnf$Sample) %>% prop.table()
  # sizes
  df.size <- aggregate(data.frame(area = assays(rnf)[["Nucleus_Area"]][1,]), 
                       by = list(cell_type = rnf[[cell.type.variable]]), FUN = "median")
  df.iter <- cbind(df.prop, df.size)
  df.iter <- df.iter[,c(1,2,3,5)]
  colnames(df.iter) <- c("cell_type", "sample_id", "true_proportion", "cell_size")
  df.iter
}))

#---------------
# get cell sizes
#---------------
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
                    s.set4 = s.set.osm.area, 
                    s.set5 = s.set.osm.mrna,
                    s.set6 = c("glial" = median(df.rn[df.rn$cell_type=="glial",]$cell_size),
                               "neuron" = median(df.rn[df.rn$cell_type=="neuron",]$cell_size)))

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
                           deconvolution.algorithm = deconvolution.algorithm)$deconvolution.results@predictions.table %>%
      as.data.frame()
    prop.pred.iter$s.set.label <- s.set.name
    prop.pred.iter$s.set.values <- paste0(
      names(s.vector.pred),"=",s.vector.pred,collapse=",")
    prop.pred.iter$sample.id <- sample.id
    prop.pred.iter$k.type <- celltype.variable
    # get true proportions
    df.rn.iter <- df.rn[df.rn$sample_id==sample.id,]
    prop.pred.iter$true.glial <- df.rn.iter[df.rn.iter$cell_type=="glial",]$true_proportion
    prop.pred.iter$true.neuron <- df.rn.iter[df.rn.iter$cell_type=="neuron",]$true_proportion
    prop.pred.iter
  }))
}))
df.s.k2.shared$experiment.type <- "shared.reference"

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
                           deconvolution.algorithm = deconvolution.algorithm)$deconvolution.results@predictions.table %>%
      as.data.frame()
    prop.pred.iter$s.set.label <- s.set.name
    prop.pred.iter$s.set.values <- paste0(
      names(s.vector.pred),"=",s.vector.pred,collapse=",")
    prop.pred.iter$sample.id <- sample.id
    prop.pred.iter$k.type <- celltype.variable
    # get true proportions
    df.rn.iter <- df.rn[df.rn$sample_id==sample.id,]
    prop.pred.iter$true.glial <- df.rn.iter[df.rn.iter$cell_type=="glial",]$true_proportion
    prop.pred.iter$true.neuron <- df.rn.iter[df.rn.iter$cell_type=="neuron",]$true_proportion
    prop.pred.iter
  }))
}))
df.s.k2.within$experiment.type <- "within.reference"

#---------------
# format results
#---------------
df.k2 <- rbind(df.s.k2.shared, df.s.k2.within)
# append abs.error
df.k2$abs.error.neuron <- abs(df.k2$neuron-df.k2$true.neuron)
df.k2$abs.error.glial <- abs(df.k2$glial-df.k2$true.glial)

#-----
# plot
#-----
# proportions scatterplots
ggplot(df.k2, aes(x = true.neuron, y = neuron)) + theme_bw() +
  geom_point(alpha = 0.5) + geom_abline(slope = 1, intercept = 0) + 
  facet_wrap(~experiment.type*s.set.label)

ggplot(df.k2, aes(x = true.neuron, y = neuron)) + theme_bw() +
  geom_point(alpha = 0.5) + geom_abline(slope = 1, intercept = 0) + 
  facet_wrap(~s.set.label)

# absolute errors
ggplot(df.k2, aes(x = experiment.type, y = abs.error.neuron)) + 
  geom_jitter(alpha = 0.5) + geom_boxplot(alpha = 0, color = "cyan")

ggplot(df.k2, aes(x = experiment.type, y = abs.error.neuron)) + 
  geom_jitter(alpha = 0.5) + geom_boxplot(alpha = 0, color = "cyan") +
  facet_wrap(~s.set.label) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
