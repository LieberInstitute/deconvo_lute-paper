
# helper functions
sn_eset_rescale <- function(sn.eset, s.glial = 3, s.neuron = 10){
  mexprs <- exprs(sn.eset)
  summary(mexprs[,1])
  which.glial <- sn.eset[["k2"]]=="glial"
  mexprs[,which.glial] <- mexprs[,which.glial] * s.glial
  which.neuron <- sn.eset[["k2"]]=="neuron"
  mexprs[,which.neuron] <- mexprs[,which.neuron] * s.neuron
  exprs(sn.eset) <- mexprs
  summary(exprs(sn.eset)[,1])
  sn.eset
}

prop_adj_results <- function(mae, bisque.sce){
  # prep
  sample.id <- unique(colData(mae)$sample.id)
  sce <- mae[["snrnaseq.k2.all"]]
  y <- assays(mae[["bulk.rnaseq"]])[["counts"]]
  z <- lute::get_z_from_sce(sce, "counts", "k2")
  # prep bisque
  bulk.eset <- mae[["bulk.rnaseq"]]
  bulk.eset <- se_to_eset(bulk.eset)
  bulk.eset[["sample.id"]] <- bulk.eset[["batch.id2"]]
  # prep bisque sce
  sn.eset <- sce_to_eset(bisque.sce)
  sn.eset[["sample.id"]] <- sn.eset[["Sample"]]
  s.vector.scale <- c("glial" = 3, "neuron" = 10)
  s.vector.noscale <- c("glial" = 1, "neuron" = 1)
  sn.eset.rescale <- sn_eset_rescale(sn.eset, 0.1, 100) # bisque rescale
  # experiment -- nnls
  nnls.scale <- lute(z = z, y = y, 
                     s = s.vector.scale, 
                     typemarker.algorithm = NULL)$deconvolution.results@predictions.table
  nnls.noscale <- lute(z = z, y = y, 
                       s = s.vector.noscale, 
                       typemarker.algorithm = NULL)$deconvolution.results@predictions.table
  colnames(nnls.scale) <- paste0(colnames(nnls.scale), ".nnls.scale")
  colnames(nnls.noscale) <- paste0(colnames(nnls.noscale), ".nnls.noscale")
  # experiment -- music
  music.prop.scale <- deconvolution(musicParam(y, z, s = s.vector.scale))@predictions.table
  music.prop.noscale <- deconvolution(musicParam(y, z, s = s.vector.noscale))@predictions.table
  colnames(music.prop.scale) <- paste0(colnames(music.prop.scale), ".music.scale")
  colnames(music.prop.noscale) <- paste0(colnames(music.prop.noscale), ".music.noscale")
  # experiment -- bisque
  # bisque no scale
  bisque.prop.noscale <- BisqueRNA::ReferenceBasedDecomposition(
    bulk.eset, sn.eset, markers=NULL, use.overlap=FALSE,
    subject.names = "sample.id", cell.types = "k2")$bulk.props
  bisque.prop.rescale <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, 
                                                                sn.eset.rescale, 
                                                                markers=NULL, use.overlap=FALSE,
                                                                subject.names = "sample.id",
                                                                cell.types = "k2")$bulk.props
  bisque.prop.noscale <- as.data.frame(t(bisque.prop.noscale))
  bisque.prop.rescale <- as.data.frame(t(bisque.prop.rescale))
  colnames(bisque.prop.rescale) <- paste0(colnames(bisque.prop.rescale), ".bisque.scale")
  colnames(bisque.prop.noscale) <- paste0(colnames(bisque.prop.noscale), ".bisque.noscale")
  # bind results
  df.res <- cbind(nnls.noscale, nnls.scale,
                  music.prop.noscale, music.prop.scale,
                  bisque.prop.rescale, bisque.prop.noscale)
  df.res <- as.data.frame(df.res)
  df.res$neuron.true <- prop.table(table(sce[["k2"]]))[["neuron"]]
  df.res$sample.id <- sample.id
  # get plot
  gg.pairs.plot <- ggpairs(df.res[,grepl("neuron", colnames(df.res))])
  return(list(df.res = df.res, gg.pairs.plot = gg.pairs.plot))
}

experiment_all_samples <- function(sample.id.vector, mae){
  bisque.sce <- mae[["snrnaseq.k2.all"]]
  lr <- list()
  for(sample.id in sample.id.vector){
    message("working on sample: ", sample.id)
    mae.iter <- mae[,colData(mae)$sample.id==sample.id,]
    lr[[sample.id]] <- prop_adj_results(mae.iter, bisque.sce)
  }
  
  #lr <- lapply(sample.id.vector, function(sample.id){
  #  message("working on sample: ", sample.id)
  #  mae.iter <- mae[,colData(mae)$sample.id==sample.id,]
  #  prop_adj_results(mae.iter, bisque.sce)
  #})
  names(lr) <- sample.id.vector
  return(lr)
}