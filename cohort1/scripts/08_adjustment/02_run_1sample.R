#!/usr/bin/env R

# Author: Sean Maden
#
# Performs bias-adjusted deconvolution.
#

source("./scripts/08_adjustment/00_musicParam-class.R")
source("./scripts/06_estimate/00_param.R")
source("./scripts/08_adjustment/00_param.R")


libv <- c("snow", "dplyr", "parallel", "doParallel", "lute", "BisqueRNA", "MuSiC", 
          "dplyr", "MultiAssayExperiment", "GGally")
sapply(libv, library, character.only = T)

#-----
# load
#-----
new.mae.filename <- "mae_allsamples.rda"
mae.final.filepath <- file.path("outputs", "01_mae", new.mae.filename)
mae <- get(load(mae.final.filepath))
sample.id.keep <- c("Br8325_mid", "Br2720_mid")
mae <- mae[,colData(mae)$sample.id %in% sample.id.keep,]

bulk.sample.remove <- c("2107UNHS-0293_Br2720_Mid_Nuc", "2107UNHS-0293_Br2720_Mid_Cyto")
mae[["bulk.rnaseq"]] <- mae[["bulk.rnaseq"]][,!colnames(mae[["bulk.rnaseq"]]) %in% bulk.sample.remove]

#-----------
# experiment
#-----------
sample.id.vector <- colData(mae)$sample.id
list.experiment.results <- experiment_all_samples(sample.id.vector, mae)
df.res <- as.data.frame(do.call(rbind, lapply(list.experiment.results, function(item){item$df.res})))

df.res$sample.id <- gsub("_.*", "", rownames(df.res))

rm(mae)

# save image
save.image(file = "./env/08_adjustment/01_run_script.RData")

#------------
# plot neuron
#------------

df.res.neuron <- df.res[,grepl("neuron", colnames(df.res))]
df.res.neuron$sample.id <- df.res$sample.id

#df.res.neuron <- df.res.neuron[,c(seq(2,ncol(df.res.neuron),1),1)]
ggpairs.neuron <- ggpairs(df.res.neuron, columns = rev(colnames(df.res.neuron)),
                          xlim = c(0, 1), ylim = c(0, 1))

# save
jpeg("./figures/08_adjustment/figs_pairs_neuron_2samples.jpg", 
     width = 10, height = 10, units = "in", res = 400)
ggpairs.neuron
dev.off()

ggplot(df.res.neuron, aes(x = neuron.nnls.noscale, y = neuron.nnls.scale, color = sample.id)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0) + theme_bw() +
  xlim(0, 1) + ylim(0, 1)

ggplot(df.res.neuron, aes(x = neuron.music.noscale, y = neuron.music.scale, color = sample.id)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0) + theme_bw() +
  xlim(0, 1) + ylim(0, 1)

ggplot(df.res.neuron, aes(x = neuron.bisque.noscale, y = neuron.bisque.scale, color = sample.id)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0) + theme_bw() +
  xlim(0, 1) + ylim(0, 1)


#-----
# save
#-----
save.image(file = "./env/08_adjustment/01_run_script.RData")
