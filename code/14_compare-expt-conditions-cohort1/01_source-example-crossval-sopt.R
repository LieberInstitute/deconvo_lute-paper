#!/usr/bin/env R

# Author: Sean Maden
#
# Manage a cross-validation series.
#

libv <- c("lute")
sapply(libv, library, character.only = T)

#scripts.vector <- c("00_source_example-soptimize-framework-lute.R", "00_source-deconvo-plots.R")
#sapply(scripts.vector, source)

#--------------------------
# manage experiment factors
#--------------------------
get_soptimize_data_list <- function(sce = NULL, y.eset = NULL, 
                                    y.train = NULL, y.validate = NULL,
                                    list.df.true = NULL, dfs = NULL,
                                    sample.id.vector = NULL,
                                    percent.train = 80,
                                       markers.per.type = 30,
                                       num.cells = 3000,
                                       total.genes = 90,
                                       num.samples = 5,
                                       num.types = 2,
                                       min.size = 1,
                                       max.size = 10,
                                       size.step = 1,
                                       assay.name = "counts",
                                       celltype.variable = "celltype",
                                    group.name = "sample.id", 
                                    deconvolution.algorithm = "nnls",
                                    with.parallel = TRUE,
                                    matched.sce = FALSE,
                                    seed.num = 0){
  #
  # seed.num = 0
  # set.seed(seed.num)
  # list.example.data <- get_soptimize_data_list()
  # dim(list.example.data$y.train)
  # dim(list.example.data$y.validate)
  #
  #
  
  # get sce, sample.id.vector
  if(is(sce, "NULL")){
    sce <- random_sce(num.genes = total.genes, num.cells = num.cells, num.types = num.types)
    sce[["sample.id"]] <- rep(paste0("sample", seq(num.samples)), num.cells/num.samples)
    # table(sce.example$sample.id, sce.example$celltype) # check
  }
  if(is(sample.id.vector, "NULL")){
    sample.id.vector <- unique(sce$sample.id)  
  }
  if(is(y.eset, "NULL")){
    y.all <- do.call(cbind, lapply(seq(num.samples), function(index){
      sce.example <- random_sce(num.genes = total.genes)
      ypb_from_sce(sce = sce, assay.name = assay.name, 
                   celltype.variable = celltype.variable)
    }))
    # format bulk eset
    colnames(y.all) <- paste0("sample", seq(ncol(y.all)))
    pdata <- data.frame(sample.id = colnames(y.all))
    rownames(pdata) <- colnames(y.all)
    y.eset <- ExpressionSet(assayData = as.matrix(y.all), 
                            phenoData = AnnotatedDataFrame(pdata))
    y.eset <- eset_to_se(y.eset)
    # subset cross-validation groups
    index.all <- seq(num.samples)
    index.train <- round((percent.train*num.samples)/100, 0)
    index.train <- seq(index.train)
    index.validate <- index.all[!index.all %in% index.train]
    y.train <- y.eset[,seq(index.train)]
    y.validate <- y.eset[,index.validate]
  }
  # data.frame of s scale factors by cell types
  if(is(dfs, "NULL")){
    dfs <- get_dfs(num.types, min.size, max.size, size.step)
  }
  # list of true cell type proportions, by samples
  if(is(list.df.true, "NULL")){
    list.df.true <- get_list_df_true(num.samples = num.samples, 
                                     num.types = num.types, prop.step = 1e-1)
  }
  # return
  metadata <- list(percent.train = percent.train,
                   markers.per.type = markers.per.type,
                   num.cells = num.cells,
                   total.genes = total.genes,
                   num.samples = num.samples,
                   num.types = num.types,
                   min.size = min.size,
                   max.size = max.size,
                   size.step = size.step,
                   assay.name = assay.name,
                   celltype.variable = celltype.variable,
                   group.name = group.name,
                   deconvolution.algorithm = deconvolution.algorithm,
                   with.parallel = with.parallel,
                   seed.num = seed.num,
                   matched.sce = matched.sce)
  lr <- list(sce = sce, 
             sample.id.vector = sample.id.vector,
             list.df.true = list.df.true, 
             dfs = dfs,
             y.eset = y.eset,
             y.train = y.train,
             y.validate = y.validate,
             metadata = metadata)
  return(lr)
}

get_dfs <- function(num.types = 2, min.size = 1, max.size = 10, size.step = 0.1){
  text.string <- paste0(rep("seq(min.size, max.size, size.step)", num.types), collapse = ",")
  dfs <- eval(parse(text = paste0("expand.grid(", text.string, ")")))
  colnames(dfs) <- paste0("s.type", seq(ncol(dfs)))
  return(dfs)
}

# list of true cell type proportions, by samples
get_list_df_true <- function(num.samples = 2, num.types = 2, prop.step = 1e-1){
  df.prop <- get_dfs(num.types = num.types, min.size = 0, max.size = 1, size.step = prop.step)
  df.prop <- df.prop[sample(seq(nrow(df.prop)), size = num.samples, replace = T),]
  df.prop <- df.prop/rowSums(df.prop)
  colnames(df.prop) <- paste0("type", seq(ncol(df.prop)))
  list.df.true <- lapply(seq(nrow(df.prop)), function(iter){df.prop[iter,]})
  names(list.df.true) <- paste0("sample", seq(length(list.df.true)))
  return(list.df.true)
}

#-----------------------------
# get cross-validation results
#-----------------------------
crossvalidate_train <- function(data, s.step.validate = 10){
  #
  # crossvalidate_train
  # data : experiment data
  # s.step.validate: s.step value for dfs.validate
  #
  
  # parse train
  message("working on test type: training")
  df.res.all <- tryCatch({
    df.res <- multigroup_bias_matched(sample.id.vector = data$sample.id.vector, 
                                      deconvolution.algorithm = data$metadata$deconvolution.algorithm,
                                      list.df.true = data$list.df.true, 
                                      y.unadj = data$y.train, 
                                      dfs = data$dfs, 
                                      sce = data$sce, 
                                      assay.name = data$metadata$assay.name,
                                      celltype.variable = data$metadata$celltype.variable,
                                      y.group.name = data$metadata$group.name,
                                      sce.group.name = data$metadata$group.name,
                                      with.parallel = data$metadata$with.parallel,
                                      matched.sce = data$metadata$matched.sce)
    df.res$test.type <- "training"
    df.res
  })
  df.res.train <- as.data.frame(df.res.all)
  
  message("updating s search space...")
  which.error.column <- which(grepl('^error\\..*', colnames(df.res.all)))[1]
  error.vector <- df.res.all[,which.error.column]
  min.decile.error <- as.numeric(quantile(error.vector)[[2]])
  error.dec.filter <- error.vector <= min.decile.error
  df.res.min <- df.res.all[error.dec.filter,]
  min.s <- min(unlist(df.res.min[,grepl("^s\\..*", colnames(df.res.min))]))
  max.s <- max(unlist(df.res.min[,grepl("^s\\..*", colnames(df.res.min))]))
  unique.types <- unique(data$sce[[data$metadata$celltype.variable]])
  num.types <- length(unique.types)
  dfs.validate <- get_dfs(num.types, 
                          min.size = min.s, 
                          max.size = max.s, 
                          size.step = (max.s-min.s)/s.step.validate)
  colnames(dfs.validate) <- paste0("s.", unique.types)
  
  lr <- list(df.res.train = df.res.train, dfs.validate = dfs.validate)
  return(lr)
}

crossvalidate_validate <- function(data, dfs.validate){
  #
  # crossvalidate_validate
  # data : experiment data object
  # dfs.validate : dfs.validate returned by crossvalidate_train()
  #
  #
  
  # parse validate
  message("working on test type: validation")
  df.res.all <- tryCatch({
    df.res <- multigroup_bias_matched(sample.id.vector = data$sample.id.vector, 
                                      deconvolution.algorithm = data$metadata$deconvolution.algorithm,
                                      list.df.true = data$list.df.true, 
                                      y.unadj = data$y.validate, 
                                      dfs = dfs.validate, 
                                      sce = data$sce, 
                                      assay.name = data$metadata$assay.name,
                                      celltype.variable = data$metadata$celltype.variable,
                                      y.group.name = data$metadata$group.name,
                                      sce.group.name = data$metadata$group.name,
                                      with.parallel = data$metadata$with.parallel,
                                      matched.sce = data$metadata$matched.sce)
    df.res$test.type <- "validation"
    df.res
  })
  df.res.validate <- as.data.frame(df.res.all)
  lr <- list(df.res.validate = df.res.validate)
  return(lr)
}

crossvalidate_soptimization <- function(list.soptimize.data, 
                                        s.step.validate = 10,
                                        facet.variable = NULL, 
                                        update.stat.summaries = TRUE,
                                        draw.min.err.line = TRUE,
                                        plot.results = FALSE){
  #
  #
  # crossvalidate_soptimization
  #
  #
  # list.soptimize.data <- get_soptimize_data_list(num.types = 3, size.step = 2)
  # crossvalidate.results <- crossvalidate_soptimization(list.soptimize.data)
  # names(crossvalidate.results)
  # crossvalidate.results$plot.list.results$`type1;type2`$heatmaps$heatmap1
  # crossvalidate.results$plot.list.results$`type1;type2`$heatmaps$heatmap6
  #
  
  data <- list.soptimize.data
  train.result <- crossvalidate_train(data, s.step.validate)
  validate.result <- crossvalidate_validate(data, train.result[["dfs.validate"]])
  # return list
  lr <- list(df.res.train = train.result[["df.res.train"]],
             df.res.validate = validate.result[["df.res.validate"]])
  if(plot.results){
    # get results plots
    lr[["plot.res.train"]] <- tryCatch({
      deconvo_plots_list(df.res.train, facet.variable, 
                         update.stat.summaries, draw.min.err.line)
    })
    lr[["plot.res.validate"]] <- tryCatch({
      deconvo_plots_list(df.res.validate, facet.variable, 
                         update.stat.summaries, draw.min.err.line)
      })
  }
  return(lr)
}

#-------------------------------
# get k match experiment results
#-------------------------------
kmatch_experiment <- function(k.variable.name = "k2",
                              sce = sce, 
                              dfs.train = NULL,
                              sample.id.vector = sample.id.vector, 
                              list.df.true = list.df.true, 
                              y.eset = y.unadj, y.train = y.train, 
                              y.validate = y.validate, 
                              assay.name = assay.name, 
                              group.name = group.name,
                              plot.option = FALSE){
  k.cell.type.vector <- unique(sce[[k.variable.name]])
  num.types <- length(k.cell.type.vector)
  if(is(dfs.train, "NULL")){
    # define dfs.train
    dfs.train <- get_dfs(num.types, min.size = 1, max.size = 5, size.step = 1)
  }
  colnames(dfs.train) <- paste0("s.", k.cell.type.vector)
  if(is(list.df.true, "NULL")){
    y.sample.id.vector <- unique(y.unadj[[y.group.variable.name]])
    list.df.true <- df.true.list(df.rn, y.sample.id.vector, 
                                 k.variable.name, k.cell.type.vector)
    names(list.df.true) <- y.sample.id.vector
  }
  # Get experiment data
  list.expt <- get_soptimize_data_list(sce = sce, 
                                       sample.id.vector = sample.id.vector, 
                                       list.df.true = list.df.true, dfs = dfs.train,
                                       y.eset = y.unadj, y.train = y.train, 
                                       y.validate = y.validate, 
                                       assay.name = assay.name, 
                                       celltype.variable = k.variable.name, 
                                       group.name = group.name,
                                       matched.sce = FALSE)
  # get results
  message("getting matched results...")
  crossval.unmatched.result <- crossvalidate_soptimization(list.expt, plot.results = plot.option)
  message("getting unmatched results...")
  list.expt$metadata$matched.sce <- TRUE
  crossval.matched.result <- crossvalidate_soptimization(list.expt, plot.results = plot.option)
  # return
  lr <- list(list.expt = list.expt,
             crossval.unmatched.result = crossval.unmatched.result,
             crossval.matched.result = crossval.matched.result)
  return(lr)
}

kmatch_experiment_plots <- function(kmatch.expt.result){
  # inspect results
  df.matched <- kmatch.expt.result$crossval.matched.result$df.res.validate
  df.unmatched <- kmatch.expt.result$crossval.unmatched.result$df.res.validate
  
  # summarize
  summary.matched.error <- summary(df.matched$error.neuron.true.pred)
  summary.unmatched.error <- summary(df.unmatched$error.neuron.true.pred)
  
  # plots -- tall
  dfp.matched <- kmatch.expt.result$crossval.matched.result$plot.list.results[[1]]$dfp
  dfp.unmatched <- kmatch.expt.result$crossval.unmatched.result$plot.list.results[[1]]$dfp
  dfp.matched$group <- "matched"
  dfp.unmatched$group <- "unmatched"
  dfp <- rbind(dfp.matched, dfp.unmatched)
  filter.dfp <- dfp$all.highlight.categories=="min.dec"
  dfp <- dfp[filter.dfp,]
  # violin plots
  vp1 <- ggplot(dfp, aes(x = sample.id, y = error.value)) + 
    geom_violin(draw_quantiles = 0.5) + facet_wrap(~group, nrow = 2) + 
    theme(axis.text.x = element_text(angle = 45,hjust = 1))
  # jitterbox plots
  jb1 <- ggplot(dfp, aes(x = sample.id, y = error.value)) + 
    geom_jitter(alpha = 0.5) + geom_boxplot(alpha = 0, color = "cyan") + 
    facet_wrap(~group, nrow = 2) + 
    theme(axis.text.x = element_text(angle = 45,hjust = 1))
  
  # plots -- wide
  dfp <- data.frame(error.matched = dfp.matched$error.type1,
                    error.unmatched = dfp.unmatched$error.type1)
  dfp$sample.id <- dfp.matched$sample.id
  
  pt1 <- ggplot(dfp, aes(x = error.matched, y = error.unmatched)) + 
    geom_abline(slope = 1, intercept = 0) + geom_point() + facet_wrap("sample.id")
  
  # filter min.dec.err.match
  dfp$min.dec.err.match <- dfp.matched$minimum.decile.error
  dfpf <- dfp[dfp$min.dec.err.match==T,]
  pt2 <- ggplot(dfpf, aes(x = error.matched, y = error.unmatched)) + 
    geom_abline(slope = 1, intercept = 0) + geom_point() + 
    facet_wrap(~sample.id)
  
  # filter min.dec.err.unmatch
  dfp$min.dec.err.unmatch <- dfp.unmatched$minimum.decile.error
  dfpf <- dfp[dfp$min.dec.err.unmatch==T,]
  pt3 <- ggplot(dfpf, aes(x = error.matched, y = error.unmatched)) + 
    geom_abline(slope = 1, intercept = 0) + geom_point() + 
    facet_wrap(~sample.id)
  
  # composite
  mindecerr.matched <- dfp.matched$minimum.decile.error
  mindecerr.unmatched <- dfp.unmatched$minimum.decile.error
  pt4 <- dfp$category <- ifelse(mindecerr.matched & mindecerr.unmatched, "both",
                         ifelse(mindecerr.matched, "matched",
                                ifelse(mindecerr.unmatched, "unmatch", "high")))
  
  pt5 <- ggplot(dfp, aes(x = error.matched, y = error.unmatched, 
                  group = category, color = category)) + 
    geom_abline(slope = 1, intercept = 0) + geom_point() + 
    facet_wrap(~sample.id)
  
  pt6 <- ggplot(dfp, aes(x = error.matched, y = error.unmatched, 
                  group = category, color = category)) + 
    geom_abline(slope = 1, intercept = 0) + geom_point() + 
    facet_wrap(~category)
  
  pt7 <- ggplot(dfp[dfp$category=="both",], aes(x = error.matched, y = error.unmatched, 
                                         group = category, color = category)) + 
    geom_abline(slope = 1, intercept = 0) + geom_point() + 
    facet_wrap(~sample.id)
  
  lr <- list(summary.matched.error = summary.matched.error,
             summary.unmatched.error = summary.unmatched.error,
             vp1 = vp1, jb1 = jb1, pt1 = pt1, pt2 = pt2,
             pt3 = pt3, pt4 = pt4, pt5 = pt5, pt6 = pt6,
             pt7 = pt7)
  return(lr)
}

#-------------------------
# s optimization functions
#-------------------------
# parallel_bias_matched
# get bias computations in parallel (THIS SCRIPT, AND A FEW OTHERS)
parallel_bias_matched <- function(sce.iter, y.iter, dfs, 
                                  celltype.variable, 
                                  deconvolution.algorithm = "nnls",
                                  df.true = NULL, 
                                  assay.name = "counts",
                                  with.parallel = TRUE){
  if(with.parallel){
    # begin parallel
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
  }
  # parse df.true options
  if(is(df.true, "NULL")){
    df.true <- sce[[celltype.variable]] %>% 
      table() %>% prop.table() %>% as.data.frame()
    rownames(df.true) <- df.true[,1]
  }
  # get full run
  if(is(y.iter, "RangedSummarizedExperiment")|is(y.iter, "SummarizedExperiment")){
    y.iter <- assays(y.iter)[[assay.name]] %>% as.matrix()
  }
  # get first s vector
  s.vector <- unlist(as.vector(dfs[1, grepl("^s\\..*", colnames(dfs))]))
  names(s.vector) <- gsub("s.", "", names(s.vector))
  # get results data.frame
  df.res <- do.call(rbind, 
                    mclapply(seq(nrow(dfs)), 
                             function(i){
                               s.vector <- unlist(
                                 as.vector(
                                   dfs[i,grepl("^s\\..*", colnames(dfs))]))
                               names(s.vector) <- gsub("s.", "", names(s.vector))
                               suppressMessages(
                                 dfi <- lute(sce.iter, y = y.iter, 
                                             celltype.variable = celltype.variable, 
                                             s = s.vector, assay.name = assay.name, 
                                             deconvolution.algorithm = deconvolution.algorithm,
                                             typemarker.algorithm = NULL)$deconvolution.results@predictions.table
                               )
                               dfi$y.sample.label <- colnames(y.iter)
                               dfi <- cbind(dfi, dfs[i,]) # append all dfs data
                               return(dfi)
                             }))
  # format results
  colnames(df.res)[1:length(s.vector)] <- paste0(
    colnames(df.res)[1:length(s.vector)], ".pred.",deconvolution.algorithm)
  # append new columns by type
  for(type.iter in names(s.vector)){
    # append true
    df.res[,ncol(df.res)+1] <- as.numeric(df.true[type.iter,1])
    colnames(df.res)[ncol(df.res)] <- paste0(type.iter, ".true")
    # append bias
    which.pred <- grepl(paste0("^", type.iter, "\\.pred\\.."), colnames(df.res))
    df.res[,ncol(df.res)+1] <- df.res[,ncol(df.res)] - df.res[,which.pred]
    colnames(df.res)[ncol(df.res)] <- paste0("bias.", type.iter, ".true.pred")
    # append error
    df.res[,ncol(df.res)+1] <- abs(df.res[,ncol(df.res)])
    colnames(df.res)[ncol(df.res)] <- paste0("error.", type.iter, ".true.pred")
  }
  if(with.parallel){
    # make sequential again (i.e. cancels parallel)
    registerDoSEQ()
    stopCluster(cl)
  }
  return(df.res)
}

# multigroup_bias_matched
# wraps parallel_bias_matched for multiple groups, uses df.true.list
multigroup_bias_matched <- function(sample.id.vector, list.df.true, y.unadj, dfs, sce,
                                    deconvolution.algorithm = "nnls", matched.sce = FALSE,
                                    y.group.name = "sample.id", sce.group.name = "sample.id",
                                    celltype.variable = "celltype", assay.name = "counts",
                                    with.parallel = TRUE, max.markers = 1000, markers.per.type = 20){
  # parse sce
  if(nrow(sce) > max.markers){
    message("sce exceeds max markers")
    message("getting ",markers.per.type," markers per type...")
    marker.vector <- lute(sce, markers.per.type = markers.per.type,
                          celltype.variable = celltype.variable,
                          assay.name = assay.name,
                          deconvolution.algorithm = NULL)$typemarker.results
    sce <- sce[rownames(sce) %in% marker.vector,]
  }
  # filter sample.id.vector on y bulk sample ids
  filter.id.vector <- sample.id.vector %in% unique(y.unadj[[y.group.name]])
  sample.id.vector <- sample.id.vector[filter.id.vector]
  df.res <- do.call(rbind, lapply(sample.id.vector, function(sample.id){
    message(sample.id)
    # parse sce.iter
    if(matched.sce){
      sce.iter <- sce[,sce[[sce.group.name]]==sample.id]
    } else{
      sce.iter <- sce
    }
    if(nrow(sce.iter) > max.markers){
      message("sce.iter exceeds max markers")
      message("getting ",markers.per.type," markers per type...")
      marker.vector <- lute(sce.iter, markers.per.type = markers.per.type,
                            celltype.variable = celltype.variable,
                            assay.name = assay.name,
                            deconvolution.algorithm = NULL)$typemarker.results
      sce.iter <- sce.iter[rownames(sce.iter) %in% marker.vector,]
    }
    if(ncol(sce.iter)==0){return()}
    y.iter <- y.unadj[,colData(y.unadj)[,y.group.name]==sample.id]
    if(ncol(y.iter)==0){return()}
    if(sample.id %in% names(list.df.true)){
      df.true.iter <- list.df.true[[sample.id]]
    } 
    if(is(df.true.iter, "NULL")){
    } else{
      df.true.iter <- df.true.iter %>% t() %>% as.data.frame()
    }
    tryCatch({
      df.res.iter <- parallel_bias_matched(sce.iter, y.iter, dfs, 
                                           df.true = df.true.iter, 
                                           deconvolution.algorithm = deconvolution.algorithm,
                                           celltype.variable = celltype.variable, 
                                           assay.name = assay.name,
                                           with.parallel = with.parallel)
      df.res.iter$sample.id <- sample.id
    }, finally = return(df.res.iter))
  }))
  return(df.res)
}

