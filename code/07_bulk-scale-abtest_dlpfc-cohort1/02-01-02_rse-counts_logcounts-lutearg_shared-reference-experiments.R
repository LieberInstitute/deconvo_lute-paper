#---------------------------------------------------
# k2 experiment -- same reference across experiments
#---------------------------------------------------
# run k2 experiment
assay.name <- "logcounts"
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
    filter.y.sample <- colData(rse.counts)[,"batch.id2"]==sample.id
    filter.y.marker <- rownames(rse.counts) %in% rownames(sce.iter)
    rse.iter <- rse.counts[filter.y.marker, filter.y.sample]
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
    prop.pred.iter$assay.name.lutearg <- assay.name
    # get true proportions
    df.rn.iter <- df.rn[df.rn$sample_id==sample.id,]
    prop.pred.iter$true.glial <- df.rn.iter[df.rn.iter$cell_type=="glial",]$true_proportion
    prop.pred.iter$true.neuron <- df.rn.iter[df.rn.iter$cell_type=="neuron",]$true_proportion
    prop.pred.iter
  }))
}))
df.s.k2.shared$experiment.type <- "shared.reference"
df.s.k2.shared$bulk.scale.type <- "counts"
df.s.k2.shared.counts.logcounts <- df.s.k2.shared
rm(df.s.k2.shared)