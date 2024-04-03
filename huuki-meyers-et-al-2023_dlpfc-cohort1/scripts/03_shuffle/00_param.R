
libv <- c("snow", "dplyr", "parallel", "doParallel", "lute", "dplyr")
sapply(libv, library, character.only = T)

#------
# setup
#------
# helper functions
# get series of s cell size factors (THIS SCRIPT, AND A FEW OTHERS)
dfs.series <- function(s.glial.series = seq(1, 10, 0.2)){
  s.neuron.series <- rev(s.glial.series)
  dfs.series <- do.call(rbind, lapply(seq(length(s.glial.series)), function(index1){
    do.call(rbind, lapply(seq(length(s.neuron.series)), function(index2){
      c("glial" = s.glial.series[index1], "neuron" = s.neuron.series[index2])
    }))
  })) %>% as.data.frame()
  #plot(dfs.series$glial, dfs.series$neuron)
  return(dfs.series)
}

# parallel_bias_matched
# get bias computations in parallel (THIS SCRIPT, AND A FEW OTHERS)
parallel_bias_matched <- function(sce.iter, y.iter, dfs, df.true = NULL, 
                                  celltype.variable = "k2", assay.name = "counts", 
                                  s.vector.ypb = c("glial" = 3, "neuron" = 10)){
  # begin parallel
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  # get full run
  if(is(y.iter, "RangedSummarizedExperiment")|is(y.iter, "SummarizedExperiment")){
    y.iter <- assays(y.iter)[[assay.name]] %>% as.matrix()
  }
  df.res <- do.call(rbind, 
                    mclapply(seq(nrow(dfs)), 
                             function(i){
                               s.vector <- c("glial" = dfs$s.glial[i], "neuron" = dfs$s.neuron[i])
                               suppressMessages(
                                 dfi <- lute(singleCellExperiment = sce.iter, 
                                             bulkExpression = y.iter, 
                                             cellTypeVariable = celltype.variable, 
                                             cellScaleFactors = s.vector, 
                                             assayName = assay.name, 
                                             typemarkerAlgorithm = NULL)$deconvolution.results@predictions.table
                               )
                               dfi$sample.label <- colnames(y.iter)
                               dfi <- cbind(dfi, dfs[i,]) # append all dfs data
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

# multigroup_bias_matched
# wraps parallel_bias_matched for multiple groups, uses df.true.list
multigroup_bias_matched <- function(sample.id.vector, df.true.list, y.unadj, dfs, sce, 
                                    y.group.name = "batch.id2", sce.group.name = "Sample",
                                    celltype.variable = "k2", assay.name = "counts"){
  df.res <- do.call(rbind, lapply(sample.id.vector, function(sample.id){
    message(sample.id)
    sce.iter <- sce[,sce[[sce.group.name]]==sample.id]
    y.iter <- y.unadj[,colData(y.unadj)[,y.group.name]==sample.id]
    df.true.iter <- df.true.list[[sample.id]] %>% t() %>% as.data.frame()
    #colnames(df.true.iter.transpose) <- colnames(df.true.list[[sample.id]])
    df.res.iter <- parallel_bias_matched(sce.iter, y.iter, dfs, df.true.iter, 
                                         celltype.variable, assay.name)
    df.res.iter$sample.id <- sample.id
    return(df.res.iter)
  }))
  return(df.res)
}


# get true cell proportions info
rnascope_cell_info <- function(df.rn, sample.id = NULL, k.type = NULL, cell.types = NULL){
  #
  #
  # get cell info data from summary df
  #
  # example:
  # df.rn <- mae[["df.cellstat.rnascope"]]
  # ct.info <- rnascope_cell_info(df.rn, "k2", c("glial", "neuron"))
  #
  filter.vector <- colnames(df.rn)
  filter.condition.vector <- rep(TRUE, length(filter.vector))
  ct.string <- paste0(cell.types, collapse = "|")
  if(!is(ct.string, "NULL")){
    filter.condition.vector <- filter.condition.vector & grepl(ct.string, filter.vector)
  }
  if(!is(sample.id, "NULL")){
    filter.condition.vector <- filter.condition.vector & grepl(sample.id, filter.vector)
  }
  if(!is(k.type, "NULL")){
    filter.condition.vector <- filter.condition.vector & grepl(k.type, filter.vector)
  }
  df.ct <- df.rn[,filter.condition.vector]
  return(df.ct)
}

df.true.list <- function(df.rn, sample.id.vector, k.type, cell.types, info = "true_proportion"){
  #
  #
  # example:
  # sample.id.vector <- unique(y.unadj$batch.id2)
  # list.dftrue <- df.true.list(df.rn, sample.id.vector, "k2", c("glial", "neuron"))
  # names(list.dftrue) <- sample.id.vector
  #
  # NOTE: does not check cell type label order .. expect 1. glial, 2. neuron !!!
  #
  
  list.dfinfo <- lapply(sample.id.vector, function(sample.id){
    rnascope_cell_info(df.rn, sample.id, k.type, cell.types)
  })
  list.dfinfo <- lapply(list.dfinfo, function(dfinfo){
    dfinfo <- as.data.frame(dfinfo["true_proportion",,drop=F])
    dfinfo <- dfinfo[,seq(length(cell.types))]
    colnames(dfinfo) <- cell.types
    return(dfinfo)
  })
  names(list.dfinfo) <- sample.id.vector
  return(list.dfinfo)
}

dfs_byvariable <- function(df.min, variable.name.vector){
  #
  # gets the s cell scale factor summaries (medians) by variables and labels
  #
  
  dfs.new <- do.call(rbind, lapply(variable.name.vector, function(variable.name){
    unique.labels <- unique(df.min[,variable.name])  
    
    dfs.iter <- do.call(rbind, lapply(unique.labels, function(label.iter){
      df.iter <- df.min[df.min[,variable.name]==label.iter,]
      matrix(c(median(df.iter[,"s.glial"]), 
               median(df.iter[,"s.neuron"]), 
               label.iter), nrow = 1)
    })) %>% as.data.frame()
    colnames(dfs.iter) <- c("s.glial", "s.neuron", "label")
    dfs.iter$variable.name <- variable.name
    return(dfs.iter)
  })) %>% as.data.frame()
  return(dfs.new)
}

condition_comparison_boxplots <- function(variable.name, variable.label, df.res.samples){
  #
  # gets box plot comparisons by condition and label in dfs, from a df.res object.
  #
  
  #variable.name <- "library.preparation"
  #variable.label <- "polyA"
  ggtitle.string <- paste0(variable.name,"==",variable.label)
  
  # filter df.res
  filter.df.res <- df.res.samples$dfs.condition.label==variable.label
  filter.df.res <- filter.df.res & df.res.samples$dfs.condition.variable.name==variable.name
  df.res <- df.res.samples[filter.df.res,]
  
  # get dfp by filter type
  # build condition filter
  filter.condition <- df.res[,colnames(df.res) == variable.name] == variable.label
  # get condition and !condition
  dfp1 <- df.res[filter.condition,]
  dfp2 <- df.res[!filter.condition,]
  dfp1$type <- "condition"
  dfp2$type <- "other"
  dfp <- rbind(dfp1, dfp2)
  
  # new plot
  #ggplot(dfp, aes(x = type, y = error.neuron)) + 
  #  geom_violin(draw_quantiles = 0.5) + ggtitle(ggtitle.string)
  ggplot(dfp, aes(x = type, y = error.neuron)) + 
    geom_jitter(alpha = 0.5) + geom_boxplot(alpha = 0, color = "cyan") +
    ggtitle(ggtitle.string)
}

get_dfcond <- function(df.res.filtered){
  #
  # gets condition-matched and -unmatched validation error summaries
  #
  
  # condition filter
  filter.cond <- gsub(";.*", "", condition.vector.train) %in% colnames(df.res.filtered)
  condition.vector.train <- condition.vector.train[filter.cond]
  # get condition summaries
  df.cond <- do.call(rbind, lapply(condition.vector.train, function(condition.iter){
    message(condition.iter)
    # get validation condition status from train s factor labels
    dfres.filter.train <- df.res.filtered[df.res.filtered$s.train.condition==condition.iter,]
    train.varname <- unique(dfres.filter.train$s.train.variable.name)
    # get validation condition filter
    train.varlabel <- unique(dfres.filter.train$s.train.variable.label)
    validation.match.filter <- dfres.filter.train[,train.varname]==train.varlabel
    # get return df
    data.frame(condition = condition.iter,
               median.val.matched.error = 
                 median(dfres.filter.train$error.neuron[validation.match.filter]),
               median.val.unmatched.error = 
                 median(dfres.filter.train$error.neuron[!validation.match.filter]))
  }))
  df.cond
}