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

# plot dfp.wide
ggplot(list.dfp[["dfp.wide"]], aes(x = noscale, y = scale, color = sample.id)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0) + facet_wrap(~algorithm) +
  xlab("No scale") + ylab("Scaled") + xlim(0, 1) + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# plot dfp.tall
ggplot(list.dfp[["dfp.tall"]], aes(x = true, y = data, color = sample.id)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0) + facet_wrap(~scale*algorithm) +
  xlab("True") + ylab("Predicted") + xlim(0, 1) + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
