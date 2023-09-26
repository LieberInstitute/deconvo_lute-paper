#!/usr/bin/env R

# Author: Sean Maden
#
# Get S estimates from bulk validation samples, as a result of multi-step S optimization.
#

libv <- c("snow", "dplyr", "parallel", "doParallel", "lute", "dplyr", "MultiAssayExperiment")
sapply(libv, library, character.only = T)

# sets variables
script.path <- file.path("scripts", "06_estimate", "00_param.R")
source(script.path)
source("./scripts/06_estimate/00_dataset_summaries.R")
source("./scripts/06_estimate/00_deconvo_plots.R")

#----------
# load data
#----------
# load mae (SEE CODE 01 OUTPUTS)
new.mae.filename <- "mae_analysis_append.rda"
mae.final.filepath <- file.path("outputs", "01_mae", new.mae.filename)
mae <- get(load(mae.final.filepath))

# true proportions from sce data
sce.all <- mae[["snrnaseq.k2.all"]]
list.df.true <- metadata(sce.all)[["list.df.true.k2"]]

#-------------------------------
# subset validation and training
#-------------------------------
samples.index <- 1
list.sample.cv <- get(load("./outputs/00_preprocess/list_snrnaseq_sampleid.rda"))
train.sample.id <- list.sample.cv[["train"]]
validate.sample.id <- list.sample.cv[["validation"]]
cd.mae <- colData(mae)
filter.mae.train <- cd.mae$sample.id %in% train.sample.id
filter.mae.validate <- cd.mae$sample.id %in% validate.sample.id
mae.train <- mae[,which(filter.mae.train)[samples.index],]
mae.validate <- mae[,which(filter.mae.validate)[samples.index],]
rm(mae)

#-----------
# experiment
#-----------
# run experiment
seq.steps.train <- seq(5, 100, 10)
list.cv <- lapply(seq.steps.train, function(train.steps){
  get_crossvalidation_results(
    mae.train, mae.validate, num.steps.train = train.steps,
    num.steps.validate = 10)
})

# get plot data
dfp <- do.call(rbind, lapply(list.cv, function(cv.result){
  c(min(cv.result$validate.result$df.res$error.neuron), cv.result$num.steps.train)
}))
dfp <- as.data.frame(dfp)
colnames(dfp) <- c("min.error.neuron.validate", "num.steps.train")

#---------------
# save plot data
#---------------
# get plot data
dfp.wide <- do.call(rbind, lapply(list.cv, function(cv.result){
  c(cv.result$num.steps.train,
    min(cv.result$validate.result$df.res$error.neuron), 
    median(cv.result$validate.result$df.res$error.neuron), 
    mean(cv.result$validate.result$df.res$error.neuron), 
    min(cv.result$train.result$df.res$error.neuron), 
    median(cv.result$train.result$df.res$error.neuron), 
    mean(cv.result$train.result$df.res$error.neuron))
}))
dfp.wide <- as.data.frame(dfp.wide)
colnames(dfp.wide) <- c("num.steps.train", 
                        "min.validate", "median.validate", "mean.validate",
                        "min.train", "median.train", "mean.train")

# get dfp.wide
dfp.tall <- do.call(rbind, lapply(list.cv, function(cv.result){
  df.iter.val <- c(
    cv.result$num.steps.train,
    min(cv.result$validate.result$df.res$error.neuron), 
    max(cv.result$validate.result$df.res$error.neuron), 
    sd(cv.result$validate.result$df.res$error.neuron), 
    median(cv.result$validate.result$df.res$error.neuron), 
    mean(cv.result$validate.result$df.res$error.neuron)
  )
  df.iter.train <- c(
    cv.result$num.steps.train,
    min(cv.result$train.result$df.res$error.neuron), 
    max(cv.result$train.result$df.res$error.neuron), 
    sd(cv.result$train.result$df.res$error.neuron), 
    median(cv.result$train.result$df.res$error.neuron), 
    mean(cv.result$train.result$df.res$error.neuron)
  )
  df.iter <- as.data.frame(rbind(df.iter.train, df.iter.val))
  colnames(df.iter) <- c("dfs.num.steps", "min", "max", "sd", "median", "mean")
  df.iter$crossvalidation <- c("train", "validate")
  return(df.iter)
}))
dfp.tall <- as.data.frame(dfp.tall)
colnames(dfp.tall) <- c("dfs.num.steps", "min", "max", "sd", "median", "mean", "crossvalidation")
list.dfp = list(dfp.wide = dfp.wide, dfp.tall = dfp.tall)
# save
save(list.dfp, file = "./outputs/06_estimate/list_results_train_plot_limited.rda")

#-----
# save
#-----
# files

list.cv.validate <- lapply(list.cv, function(item){
  list(results.validate = item$validate.result, num.steps.train = item$num.steps.train)
})
save(list.cv.validate, file = "./outputs/06_estimate/results_cv_validate_limited.rda")

save(list.cv, file = "./outputs/06_estimate/results_cv_limited.rda")
save(dfp, file = "./outputs/06_estimate/results_dfp_limited.rda")

# env
rm(mae.final)
rm(mae.validate)
rm(mae.train)
save.image("./outputs/06_estimate/02_run_limited_script.RData")
