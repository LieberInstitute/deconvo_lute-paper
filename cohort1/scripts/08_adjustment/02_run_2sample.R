#!/usr/bin/env R

# Author: Sean Maden
#
# Performs bias-adjusted deconvolution.
#

source("./scripts/08_adjustment/00_musicParam-class.R")
source("./scripts/08_adjustment/00_sopt.R")
source("./scripts/08_adjustment/00_sopt_utilities.R")
source("./scripts/08_adjustment/00_param.R")

libv <- c("snow", "dplyr", "parallel", "doParallel", "lute", "BisqueRNA", "MuSiC", 
          "dplyr", "MultiAssayExperiment", "GGally")
sapply(libv, library, character.only = T)

# params
num.dfs.steps <- 40

#-----
# load
#-----
new.mae.filename <- "mae_allsamples.rda"
mae.final.filepath <- file.path("outputs", "01_mae", new.mae.filename)
mae <- get(load(mae.final.filepath))
sample.id.keep <- c("Br8325_mid", "Br3942_mid")
mae <- mae[,colData(mae)$sample.id %in% sample.id.keep,]

#-----------
# experiment
#-----------
sample.id.vector <- colData(mae)$sample.id
list.experiment.results <- experiment_all_samples(
  sample.id.vector, mae, dfs.steps = num.dfs.steps)
df.res <- as.data.frame(
  do.call(rbind, lapply(list.experiment.results, function(item){item$df.res})))
df.res$sample.id <- gsub("_.*", "", rownames(df.res))
list.dfp <- get_dfp_list(df.res)

# save image
rm(mae)
save.image(file = "./env/08_adjustment/02_run_2sample_script.RData")
