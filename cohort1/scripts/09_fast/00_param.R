#!/usr/bin/env R

#
# Main code to run S factor optimization.
#

libv <- c("snow", "dplyr", "parallel", "doParallel", "lute", "dplyr")
sapply(libv, library, character.only = T)

#------
# setup
#------
# helper functions
# get series of s cell size factors (THIS SCRIPT, AND A FEW OTHERS)
dfs.series <- function(s.glial.series = seq(1, 20, 1)){
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

get_dfres_plots <- function(df.res, facet.variable = NULL){
  # plot params
  # source plot functions
  list.plots.dfp <- deconvo_plots_list(df.res, facet.variable)
  return(list.plots.dfp)
}

get_sopt_results <- function(mae, dfs, label = "train"){
  # set params (SEE PROJECT NOTES)
  assay.name <- "counts"
  celltype.variable <- "k2"
  sample.id.variable <- "Sample"
  y.group.name <- 'batch.id2'
  bulk.name <- "bulk.rnaseq"
  sn.name <- "snrnaseq.k2.all"
  # sample id vector
  sample.id.vector <- unique(
    intersect(
      mae[[bulk.name]][[y.group.name]], 
      mae[[sn.name]][[sample.id.variable]]))
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
                            dfs, sce.iter, 
                            assay.name = assay.name)
  })
  df.res <- do.call(rbind, lapply(list.res, function(item){item}))
  df.res$crossvalidation <- label
  # prepare and plot results
  df.res <- dfres_postprocess(df.res)
  list.plots.dfres <- get_dfres_plots(df.res)
  return(list(df.res = df.res, list.plots = list.plots.dfres))
}

get_crossvalidation_results <- function(mae.train, mae.validate, 
                                        num.steps.train = 10,
                                        num.steps.validate = 50,
                                        s.min.train = 1,
                                        s.max.train = 40){
  # train
  s.increment.train <- (s.max.train-s.min.train)/num.steps.train
  dfs <- dfs.series(
    seq(s.min.train, s.max.train, s.increment.train))
  message("beginning training")
  list.dfres.train <- get_sopt_results(mae.train, dfs, "train")
  # get dfs from train min.error coordinates
  df.res.train <- list.dfres.train$df.res
  min.error.neuron <- min(df.res.train$error.neuron)
  message("min. error train: ", min.error.neuron)
  filter.res <- df.res.train$error.neuron == min.error.neuron
  s.train.neuron <- df.res.train[filter.res,]$s.neuron
  s.train.glial <- df.res.train[filter.res,]$s.glial
  s.vector.validate <- c(s.train.neuron, s.train.glial)
  s.validate.min <- min(s.vector.validate)
  s.validate.max <- max(s.vector.validate)
  s.validate.increment <- (s.validate.max-s.validate.min)/num.steps.validate
  s.validate.seq <- seq(s.validate.min, s.validate.max, s.validate.increment)
  dfs.validate <- dfs.series(s.validate.seq)
  # validate
  message("beginning validation")
  list.dfres.validate <- get_sopt_results(mae.validate, dfs.validate, "validate")
  return(list(train.result = list.dfres.train, 
              validate.result = list.dfres.validate,
              num.steps.train = num.steps.train))
}
