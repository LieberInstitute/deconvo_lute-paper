#!/usr/bin/env R

# Author: Sean Maden
#
# Get S estimates from bulk validation samples, as a result of multi-step S optimization.
#

libv <- c("snow", "dplyr", "parallel", "doParallel", "lute", "dplyr", "MultiAssayExperiment")
sapply(libv, library, character.only = T)

# sets variables
source("source/00_sopt.R")
source("source/00_sopt_utilities.R")
source("./source/00_dataset_summaries.R")
source("./source/00_deconvo_plots.R")

#----------
# load data
#----------
# load mae (SEE CODE 01 OUTPUTS)
new.mae.filename <- "mae_analysis_append.rda"
mae.final.filepath <- file.path("outputs", "01_mae", new.mae.filename)
mae <- get(load(mae.final.filepath))
# true proportions from sce data
sn.assay.name <- "snrnaseq.k2.all"
sce.all <- mae[[sn.assay.name]]
list.df.true <- metadata(sce.all)[["list.df.true.k2"]]

#-------------------------------
# subset validation and training
#-------------------------------
samples.index.train <- 1:5
samples.index.validate <- 1:2
list.sample.cv <- get(load("./outputs/00_preprocess/list_snrnaseq_sampleid.rda"))
train.sample.id <- list.sample.cv[["train"]]
validate.sample.id <- list.sample.cv[["validation"]]

cd.mae <- colData(mae)
filter.mae.train <- cd.mae$sample.id %in% train.sample.id
filter.mae.validate <- cd.mae$sample.id %in% validate.sample.id

mae.train <- mae[,which(filter.mae.train)[samples.index.train],]
mae.validate <- mae[,which(filter.mae.validate)[samples.index.validate],]
rm(mae)

#-----------
# experiment
#-----------
# parameters for experiment
seq.steps.train <- seq(5, 100, 10)

# run experiment
list.cv <- lapply(seq.steps.train, function(train.steps){
  list.res <- get_crossvalidation_results(
    mae.train, mae.validate, num.steps.train = train.steps)
  list.res
})
names(list.cv) <- paste0("steps:", seq.steps.train)

#-----------------
# format plot data
#-----------------
# get plot data
dfp <- do.call(rbind, lapply(list.cv, function(cv.result){
  
  df.res.val <- cv.result$validate.result$df.res
  df.res.train <- cv.result$train.result$df.res
  
  c(min(df.res.val$error.neuron),
    median(df.res.val$error.neuron),
    sd(df.res.val$error.neuron),
    
    min(df.res.train$error.neuron),
    median(df.res.train$error.neuron),
    sd(df.res.train$error.neuron),
    
    length(unique(df.res.val$sample.id)),
    length(unique(df.res.train$sample.id)),
    
    cv.result$num.steps.train)
}))
dfp <- as.data.frame(dfp)
colnames(dfp) <- c("min.val", "median.val", "sd.val",
                   "min.train", "median.train", "sd.train", 
                   "num.samples.val", "num.samples.train",
                   "steps")

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

#-----
# save
#-----

# save plot data
save(list.dfp, file = "./outputs/06_estimate/list_results_train_plot_limited.rda")


# files
list.cv.validate <- lapply(list.cv, function(item){
  list(results.validate = item$validate.result, num.steps.train = item$num.steps.train)
})
save(list.cv.validate, file = "./outputs/06_estimate/results_cv_validate_limited.rda")
save(list.cv, file = "./outputs/06_estimate/results_cv_limited.rda")
save(dfp, file = "./outputs/06_estimate/results_dfp_limited.rda")

# env
rm(mae.validate)
rm(mae.train)
save.image("./outputs/06_estimate/02_run_limited_resolution_script.RData")
