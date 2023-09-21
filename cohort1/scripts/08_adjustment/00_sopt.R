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
  list.df.true <- get_df_true_list(sce, sample.id.variable, celltype.variable)
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

# multigroup_bias_matched
# wraps parallel_bias_matched for multiple groups, uses df.true.list
multigroup_bias_matched <- function(sample.id.vector, df.true.list, y.unadj, dfs, sce, 
                                    y.group.name = "batch.id2",
                                    celltype.variable = "k2", assay.name = "counts", 
                                    s.vector.ypb = c("glial" = 3, "neuron" = 10)){
  df.res <- do.call(rbind, lapply(sample.id.vector, function(sample.id){
    message(sample.id)
    y.iter <- y.unadj[,colData(y.unadj)[,y.group.name]==sample.id]
    df.true.iter <- df.true.list[[sample.id]] %>% t() %>% as.data.frame()
    #colnames(df.true.iter.transpose) <- colnames(df.true.list[[sample.id]])
    df.res.iter <- parallel_bias_matched(sce, y.iter, dfs, df.true.iter, 
                                         celltype.variable, assay.name, s.vector.ypb)
    df.res.iter$sample.id <- sample.id
    return(df.res.iter)
  }))
  return(df.res)
}

# parallel_bias_matched
# get bias computations in parallel (THIS SCRIPT, AND A FEW OTHERS)
parallel_bias_matched <- function(sce, yunadj, dfs, df.true = NULL, 
                                  celltype.variable = "k2", assay.name = "counts", 
                                  s.vector.ypb = c("glial" = 3, "neuron" = 10)){
  # begin parallel
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  # get full run
  if(is(yunadj, "RangedSummarizedExperiment")){
    yunadj <- assays(yunadj)[[1]] %>% as.matrix()
  }
  df.res <- do.call(rbind, 
                    mclapply(seq(nrow(dfs)), 
                             function(i){
                               s.vector <- c("glial" = dfs$glial[i], "neuron" = dfs$neuron[i])
                               suppressMessages(
                                 dfi <-lute(sce, y = yunadj, celltype.variable = celltype.variable, s = s.vector,
                                            typemarker.algorithm = NULL)$deconvolution.results@predictions.table
                               )
                               dfi$sample.label <- colnames(yunadj)
                               dfi$s.glial <- s.vector["glial"]
                               dfi$s.neuron <- s.vector["neuron"]
                               return(dfi)
                             }))
  colnames(df.res)[1:2] <- paste0(colnames(df.res)[1:2], ".pred.nnls")
  if(is(df.true, "NULL")){
    df.true <- sce[[celltype.variable]] %>% table() %>% prop.table() %>% as.data.frame()
    rownames(df.true) <- df.true[,1]
  }
  df.res$glial.true <- as.numeric(df.true["glial",1])
  df.res$neuron.true <- as.numeric(df.true["neuron",1])
  df.res$bias.glial.true.pred <- df.res$glial.true - df.res$glial.pred.nnls
  df.res$bias.neuron.true.pred <- df.res$neuron.true - df.res$neuron.pred.nnls
  # make sequential again (i.e. cancels parallel)
  registerDoSEQ()
  return(df.res)
}

# postprocess results in df.res
# appends metrics and statistics to the results
dfres_postprocess <- function(df.res){
  df.res.samples <- df.res
  
  # append data transformations
  # this is the chunk that sets more operants in `df.res`
  df.res.samples$s.fraction.neuron.glial <- df.res.samples$s.neuron/df.res.samples$s.glial
  df.res.samples$log.s.fraction <- log(df.res.samples$s.fraction.neuron.glial)
  df.res.samples$error.neuron <- abs(df.res.samples$bias.neuron.true.pred)
  df.res.samples$error.glial <- abs(df.res.samples$bias.glial.true.pred)
  df.res.samples$minimum.error <- df.res.samples$error.neuron==min(df.res.samples$error.neuron)
  df.res.samples$maximum.error <- df.res.samples$error.neuron==max(df.res.samples$error.neuron)
  deciles.error.neuron <- quantile(df.res.samples$error.neuron, seq(0, 1, 0.1))
  df.res.samples$minimum.decile.error <- df.res.samples$error.neuron <= deciles.error.neuron[2]
  df.res.samples$maximum.decile.error <- df.res.samples$error.neuron >= deciles.error.neuron[9]
  df.res.samples$error.neuron <- df.res.samples$bias.neuron.true.pred %>% abs()
  df.res <- df.res.samples
  # do postprocessing on the dataset
  # rename abs bias colname to error
  colnames.filter <- colnames(df.res)=="abs.bias.neuron"
  colnames(df.res)[colnames.filter] <- "error.neuron"
  # add group labels
  df.res$glial.group.label <- "glial"
  df.res$neuron.group.label <- "neuron"
  # append highlight values
  df.res$minimum.error <- df.res$error.neuron==min(df.res$error.neuron)
  df.res$maximum.error <- df.res$error.neuron==max(df.res$error.neuron)
  # assign deciles
  quant <- quantile(df.res$error.neuron, seq(0,1,0.1))
  df.res$lowest.decile.error <- df.res$error.neuron <= quant[2]
  df.res$highest.decile.error <- df.res$error.neuron>=quant[10]
  
  return(df.res)
}
