#!/usr/bin/env R

# Author: Sean Maden
#
# Run sucessive iterations of S optimization. Use the real bulk data.
#

libv <- c("MultiAssayExperiment", "ggplot2", "gridExtra")
sapply(libv, library, character.only = T)
source("./scripts/09_fast/00_param.R")
source("./source/00_dataset_summaries.R")
source("./source/00_deconvo_plots.R")
# source("./source/00_sopt.R")

# load
mae <- get(load("./outputs/01_mae/mae_analysis_append.rda"))
mae <- mae[,colData(mae)$sample.id=="Br8492_mid",]

# helper functions
dfs_iter <- function(param.iter, s.opt.iter){
  s.diff <- as.numeric(param.iter["s.diff"])
  s.step <- as.numeric(param.iter["s.step"])
  ds.iter.new <- data.frame(glial = seq(
    s.opt.iter["glial"]-s.diff,
    s.opt.iter["glial"]+s.diff,
    s.step
  ))
  ds.iter.new$neuron <- seq(
    s.opt.iter["neuron"]-s.diff,
    s.opt.iter["neuron"]+s.diff,
    s.step
  )
  return(ds.iter.new)
}

bias_res_iter <- function(sample.id.vector, df.true.list, 
                          y.unadj, dfs, sce, y.group.name,
                          celltype.variable, assay.name){
  t1 <- Sys.time()
  df.res <- multigroup_bias_matched(sample.id.vector = sample.id.vector, 
                                    df.true.list = df.true.list, 
                                    y.unadj = y.unadj, dfs = dfs, 
                                    sce = sce, y.group.name = y.group.name,
                                    celltype.variable = celltype.variable, 
                                    assay.name = assay.name)
  t2 <- Sys.time()
  time.elapsed <- t2 - t1
  df.res <- dfres_postprocess(df.res)
  filter.error <- df.res$error.neuron==min(df.res$error.neuron)
  s.opt.iter <- c(df.res[filter.error,]$s.glial,
                  df.res[filter.error,]$s.neuron)
  names(s.opt.iter) <- c("glial", "neuron")
  list.plots.res <- get_dfres_plots(df.res)
  list.res <- list(df.res = df.res, s.optimal = s.opt.iter, 
    time.elapsed = time.elapsed, list.plots.res = list.plots.res)
  return(list.res)
}

run_sopt_series <- function(dfs.param, sample.id.vector, df.true.list, y.unadj, 
                            sce, y.group.name = "batch.id2",
                            celltype.variable = "k2", 
                            assay.name = "logcounts"){
  list.res <- list()
  s.opt.iter <- NA
  for(iter in seq(nrow(dfs.param))){
    message("working on iteration ", iter, "...")
    new.res.name <- paste0("iter", iter)
    message("s.opt.iter: ", s.opt.iter)
    param.iter <- dfs.param[iter,]
    if(is(param.iter[1], "NULL")|is.na(param.iter[1])){
      if(is(s.opt.iter, "NULL")){stop("Error in dfs.param.")}
      dfs.iter.new <- dfs_iter(param.iter, s.opt.iter)
      message("num tests ", nrow(dfs.iter.new))
    } else{
      param.iter <- as.numeric(param.iter)
      seq.dfs <- seq(
        param.iter[1], param.iter[2], param.iter[3])
      dfs.iter.new <- dfs.series(seq.dfs)
    }
    message("num tests : ", nrow(dfs.iter.new), "...")
    list.res[[new.res.name]] <- 
      bias_res_iter(sample.id.vector = sample.id.vector, 
                  df.true.list = df.true.list, 
                  y.unadj = y.unadj, 
                  dfs = dfs.iter.new, 
                  sce = sce, 
                  y.group.name = y.group.name,
                  celltype.variable = celltype.variable, 
                  assay.name = assay.name)
    s.opt.iter <- list.res[[new.res.name]][["s.optimal"]]
    rm(dfs.iter.new)
  }
  return(list.res)
}

plot_sopt_series <- function(list.res){
  dfp <- do.call(rbind, lapply(seq(length(list.res)), function(index){
    res.iter <- list.res[[index]]
    data.frame(time.elapsed = as.numeric(res.iter$time.elapsed),
               min.error.neuron = min(res.iter$df.res$error.neuron),
               num.tests = nrow(res.iter$df.res),
               min.s.neuron = min(res.iter$df.res$s.neuron),
               min.s.glial = min(res.iter$df.res$s.glial),
               max.s.neuron = max(res.iter$df.res$s.neuron),
               max.s.glial = max(res.iter$df.res$s.glial),
               s.opt.neuron = res.iter$s.optimal["neuron"],
               s.opt.glial = res.iter$s.optimal["glial"],
               run.sequence = index)
  }))
  dfp <- as.data.frame(dfp)
  dfp$run.sequence <- as.numeric(dfp$run.sequence)
  # plots
  new.plot1 <- ggplot(dfp, aes(x = run.sequence, y = min.error.neuron)) + 
    geom_point() + geom_line()
  new.plot2 <- ggplot(dfp, aes(x = run.sequence, y = time.elapsed)) + 
    geom_point() + geom_line()
  new.plot3 <- ggplot(dfp, aes(x = run.sequence, y = num.tests)) + 
    geom_point() + geom_line()
  new.plot4 <- ggplot(dfp, aes(x = run.sequence, y = min.s.neuron)) + 
    geom_point() + geom_line()
  new.plot5 <- ggplot(dfp, aes(x = run.sequence, y = min.s.glial)) + 
    geom_point() + geom_line()
  new.plot6 <- ggplot(dfp, aes(x = run.sequence, y = s.opt.neuron)) + 
    geom_point() + geom_line()
  new.plot7 <- ggplot(dfp, aes(x = run.sequence, y = s.opt.glial)) + 
    geom_point() + geom_line()
  # return
  list.plots <- list(dfp = dfp,
                     min.error.neuron = new.plot1, 
                     time.elapsed = new.plot2, 
                     num.tests = new.plot3, 
                     min.s.neuron = new.plot4, 
                     min.s.glial = new.plot5,
                     s.optimal.neuron = new.plot6,
                     s.optimal.glial = new.plot7)
  return(list.plots)
}

# new run
dfs.param <- data.frame(
  s.min = c(1, NA, NA, NA),
  s.max = c(80, NA, NA, NA),
  s.step = c(2, 1e-2, 1e-3, 1e-4),
  s.diff = c(NA, 1e-1, 1e-2, 1e-3)
)
sample.id.vector <- "Br8492_mid"
df.true.list <- metadata(mae[[1]])[["list.df.true.k2"]]
y.unadj <- mae[["bulk.rnaseq"]][,1]
sce <- mae[[1]]
list.res <- run_sopt_series(dfs.param, sample.id.vector = sample.id.vector, 
                            df.true.list = df.true.list, y.unadj = y.unadj, 
                            sce = sce)

identical(list.res$iter1$df.res,
          list.res$iter2$df.res)

list.res.plots <- plot_sopt_series(list.res)

save.image(file = "./env/09_fast/01_run_bulk_1sample_script.RData")
