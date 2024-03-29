
# helper functions
sn_eset_rescale <- function(sn.eset, s.glial = 3, s.neuron = 10){
  mexprs <- m1 <- exprs(sn.eset)
  which.glial <- sn.eset[["k2"]]=="glial"
  mexprs[,which.glial] <- mexprs[,which.glial] * s.glial
  which.neuron <- sn.eset[["k2"]]=="neuron"
  mexprs[,which.neuron] <- mexprs[,which.neuron] * s.neuron
  exprs(sn.eset) <- m2 <- mexprs
  if(identical(m1, m2)){stop("rescale failure.")}
  sn.eset
}

prop_adj_results <- function(mae, bisque.sce, bisque.bulk, s.vector.scale,
                             bulk.mae.name = "bulk.rnaseq", 
                             bulk.sample.id.variable = "batch.id2",
                             assay.name = "logcounts"){
  # prep
  sample.id <- unique(colData(mae)$sample.id)
  sce <- mae[["snrnaseq.k2.all"]]
  y.set <- bisque.bulk # mae[[bulk.mae.name]]
  # get logcounts scaled expression
  y.set <- scuttle::logNormCounts(y.set)
  sce <- scuttle::logNormCounts(sce)
  bisque.sce <- scuttle::logNormCounts(bisque.sce)
  y <- assays(y.set)[["logcounts"]][,,drop=F]
  assays(sce)[["counts"]] <- assays(sce)[["logcounts"]]
  assays(bisque.sce)[["logcounts"]] <- assays(bisque.sce)[["logcounts"]]
  bulk.eset <- se_to_eset(y.set, assay.name)
  
  # prep bisque
  sn.eset <- sce_to_eset(bisque.sce)
  exprs(sn.eset) <- exprs(sn.eset) + 1e-3
  sn.eset[["sample.id"]] <- sn.eset[["Sample"]]
  bulk.eset[["sample.id"]] <- bulk.eset[[bulk.sample.id.variable]]
  # prep z
  message("prep z")
  z <- lute::referenceFromSingleCellExperiment(sce, assay.name, "k2")
  # get s vectors
  s.vector.noscale <- c("glial" = 1, "neuron" = 1)
  # bisque rescale
  message("bisque rescale")
  sn.eset.rescale <- sn_eset_rescale(sn.eset, 
                                     s.vector.scale["glial"], 
                                     s.vector.scale["neuron"]) 
  # experiment -- nnls
  message("experiment -- nnls")
  nnls.scale <- lute(singleCellExperiment = sce, 
                     bulkExpression = y, 
                     cellTypeVariable  = "k2", 
                     cellScaleFactors  = s.vector.scale, 
                     assayName = assay.name,
                     typemarkerAlgorithm = NULL)$deconvolutionResults@predictionsTable

  nnls.noscale <- lute(singleCellExperiment = sce, 
                     bulkExpression = y, 
                     cellTypeVariable  = "k2", 
                     cellScaleFactors  = s.vector.noscale, 
                     assayName = assay.name,
                     typemarkerAlgorithm = NULL)$deconvolutionResults@predictionsTable
  
  colnames(nnls.scale) <- paste0(colnames(nnls.scale), ".nnls.scale")
  colnames(nnls.noscale) <- paste0(colnames(nnls.noscale), ".nnls.noscale")
  # experiment -- music
  message("experiment -- music")
  music.prop.scale <- deconvolution(
    musicParam(bulkExpression = y, referenceExpression = z, 
               cellScaleFactors = s.vector.scale))@predictionsTable
  music.prop.noscale <- deconvolution(
    musicParam(bulkExpression = y, referenceExpression = z, 
               cellScaleFactors = s.vector.noscale))@predictionsTable
  colnames(music.prop.scale) <- 
    paste0(colnames(music.prop.scale), ".music.scale")
  colnames(music.prop.noscale) <- 
    paste0(colnames(music.prop.noscale), ".music.noscale")
  # experiment -- bisque
  message("experiment -- bisque")
  # bisque no scale
  bisque.prop.noscale <- BisqueRNA::ReferenceBasedDecomposition(
    bulk.eset, sn.eset, markers=NULL, use.overlap=FALSE,
    subject.names = "sample.id", cell.types = "k2")$bulk.props
  bisque.prop.rescale <- BisqueRNA::ReferenceBasedDecomposition(
    bulk.eset, sn.eset.rescale, markers=NULL, use.overlap=FALSE,
    subject.names = "sample.id", cell.types = "k2")$bulk.props
  
  filter.sample.id <- bulk.eset[[bulk.sample.id.variable]]==sample.id
  bisque.prop.noscale <- as.data.frame(t(bisque.prop.noscale[,filter.sample.id]))
  bisque.prop.rescale <- as.data.frame(t(bisque.prop.rescale[,filter.sample.id]))
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
  #gg.pairs.plot <- ggpairs(df.res[,grepl("neuron", colnames(df.res))])
  return(list(df.res = df.res, s.vector.scale = s.vector.scale))
}

#
# get experiments across samples group ids
experiment_all_samples <- function(sample.id.vector, mae, 
                                   list.s.vector.scale,
                                   assay.name = "logcounts",
                                   bulk.mae.name = "bulk.rnaseq", 
                                   bulk.sample.id.variable = "batch.id2"){
  bisque.sce <- mae[["snrnaseq.k2.all"]]
  bisque.bulk <- mae[[bulk.mae.name]]
  lr <- list()
  for(sample.id in sample.id.vector){
    message("working on sample: ", sample.id)
    mae.iter <- mae[,colData(mae)$sample.id==sample.id,]
    s.vector.scale <- list.s.vector.scale[[sample.id]]
    lr[[sample.id]] <- prop_adj_results(mae.iter, 
                                        bisque.sce, 
                                        bisque.bulk, 
                                        s.vector.scale = s.vector.scale,
                                        assay.name = assay.name,
                                        bulk.mae.name = bulk.mae.name, 
                                        bulk.sample.id.variable = bulk.sample.id.variable)
  }
  names(lr) <- sample.id.vector
  return(lr)
}

#
# get s optimization results table
get_sopt_results <- function(mae, dfs, label = "train",
                             assay.name = "logcounts",
                             celltype.variable = "k2",
                             sample.id.variable = "Sample",
                             y.group.name = 'batch.id2',
                             bulk.name = "bulk.rnaseq",
                             sn.name = "snrnaseq.k2.all"){
  sample.id.vector <- unique(
    intersect(
      mae[[bulk.name]][[y.group.name]], 
      mae[[sn.name]][[sample.id.variable]]))
  # get list.df.true
  sce <- mae[[sn.name]]
  # list.df.true <- get_df_true_list(sce, sample.id.variable, celltype.variable)
  list.df.true <- metadata(sce)[["list.df.true.k2"]]
  # iterate on training samples
  list.res <- lapply(sample.id.vector, function(sample.id){
    filter.mae <- colData(mae)$sample.id==sample.id
    mae.iter <- mae[,filter.mae,]
    message("Num. tests: ", nrow(dfs))
    # seq
    y.iter <- mae.iter[[bulk.name]]
    sce.iter <- mae.iter[[sn.name]]
    # results
    multigroup_bias_matched(sample.id, list.df.true, y.iter, 
                            y.group.name = y.group.name,
                            dfs, sce.iter, assay.name = assay.name)
  })
  df.res <- do.call(rbind, lapply(list.res, function(item){item}))
  df.res$crossvalidation <- label
  # prepare and plot results
  df.res <- dfres_postprocess(df.res)
  return(list(df.res = df.res))
}

# get wide and tall plot data frames
get_dfp_list <- function(df.res){
  # plot dfp.wide
  dfp.wide <- rbind(data.frame(scale = df.res$neuron.music.scale,
                               noscale = df.res$neuron.music.noscale,
                               algorithm = rep("music", nrow(df.res)),
                               sample.id = df.res$sample.id),
                    data.frame(scale = df.res$neuron.nnls.scale,
                               noscale = df.res$neuron.nnls.noscale,
                               algorithm = rep("nnls", nrow(df.res)),
                               sample.id = df.res$sample.id),
                    data.frame(scale = df.res$neuron.bisque.scale,
                               noscale = df.res$neuron.bisque.noscale,
                               algorithm = rep("bisque", nrow(df.res)),
                               sample.id = df.res$sample.id)
  )
  # plot dfp.tall
  dfp.tall <- rbind(data.frame(data = df.res$neuron.music.scale,
                               scale = rep("scale", nrow(df.res)),
                               true = df.res$neuron.true,
                               sample.id = df.res$sample.id,
                               algorithm = rep("music", nrow(df.res))),
                    data.frame(data = df.res$neuron.nnls.scale,
                               scale = rep("scale", nrow(df.res)),
                               true = df.res$neuron.true,
                               sample.id = df.res$sample.id,
                               algorithm = rep("nnls", nrow(df.res))),
                    data.frame(data = df.res$neuron.bisque.scale,
                               scale = rep("scale", nrow(df.res)),
                               true = df.res$neuron.true,
                               sample.id = df.res$sample.id,
                               algorithm = rep("bisque", nrow(df.res))),
                    data.frame(data = df.res$neuron.music.noscale,
                               scale = rep("noscale", nrow(df.res)),
                               true = df.res$neuron.true,
                               sample.id = df.res$sample.id,
                               algorithm = rep("music", nrow(df.res))),
                    data.frame(data = df.res$neuron.nnls.noscale,
                               scale = rep("noscale", nrow(df.res)),
                               true = df.res$neuron.true,
                               sample.id = df.res$sample.id,
                               algorithm = rep("nnls", nrow(df.res))),
                    data.frame(data = df.res$neuron.bisque.noscale,
                               scale = rep("noscale", nrow(df.res)),
                               true = df.res$neuron.true,
                               sample.id = df.res$sample.id,
                               algorithm = rep("bisque", nrow(df.res))))
  return(list(dfp.wide = dfp.wide, dfp.tall = dfp.tall))
}