#!/usr/bin/env R

# Author: Sean Maden
#
# Performs bias-adjusted deconvolution.
#


source("./scripts/06_estimate/00_param.R")
source("./scripts/08_adjustment/00_param.R")


libv <- c("snow", "dplyr", "parallel", "doParallel", "lute", "BisqueRNA", "MuSiC", 
          "dplyr", "MultiAssayExperiment", "GGally")
sapply(libv, library, character.only = T)

source("./scripts/08_adjustment/00_musicParam-class.R")

#-----
# load
#-----
new.mae.filename <- "mae_analysis_append.rda"
mae.final.filepath <- file.path("outputs", "01_mae", new.mae.filename)
mae <- get(load(mae.final.filepath))
sample.id.vector <- unlist(get(load("./outputs/00_preprocess/list_snrnaseq_sampleid.rda")))
mae <- mae[,colData(mae)$sample.id %in% as.character(sample.id.vector),]
nrow(colData(mae))

samples.exclude <- c("Br2743_mid", "Br6522_mid")
mae <- mae[,!colData(mae)$sample.id %in% samples.exclude,]
nrow(colData(mae))

#-----------
# experiment
#-----------
sample.id.vector <- colData(mae)$sample.id
list.experiment.results <- experiment_all_samples(sample.id.vector, mae)
df.res <- as.data.frame(
  do.call(rbind, lapply(list.experiment.results, function(item){item$df.res})))

#------------
# plot neuron
#------------
df.res.neuron <- df.res[,grepl("neuron", colnames(df.res))]
#df.res.neuron <- df.res.neuron[,c(seq(2,ncol(df.res.neuron),1),1)]
ggpairs.neuron <- ggpairs(df.res.neuron, columns = rev(colnames(df.res.neuron)))

# save
jpeg("./figures/08_adjustment/figs_pairs_neuron.jpg", 
     width = 10, height = 10, units = "in", res = 400)
ggpairs.neuron
dev.off()

#-----
# save
#-----
save.image(file = "./env/08_adjustment/01_run_script.RData")
