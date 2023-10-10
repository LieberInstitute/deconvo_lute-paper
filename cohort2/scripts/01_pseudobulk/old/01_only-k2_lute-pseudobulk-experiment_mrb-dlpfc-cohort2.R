#!/usr/bin/env R

# Author: Sean Maden
#
# Test independent pseudobulk from DLPFC cohort1 snRNAseq data.
#

source("./cohort2/scripts/01_pseudobulk/00_parameters-pseudobulk.R")
sapply(libv, library, character.only = T)

# load marker data
markers.path <- "./cohort1/outputs/00_preprocess/list-scef_markers-k2-k3-k4_ro1-dlpfc.rda"
list.sce.markers <- get(load(markers.path))

# load sce data
sce.path <- "./cohort2/data/sce-mrb_dlpfc.rda"
sce.mrb <- get(load(sce.path)) # load sce data
sce.k2 <- sce.mrb[rownames(sce.mrb) %in% rownames(list.sce.markers[["k2"]]),]

# get experiment results tables
dfp.tall <- get_ypb_experiment_series(sce.k2, sample.id.variable = "donor", 
                                      celltype.variable = "k2", assay.name = "logcounts",
                                      s.vector = c("glial" = 3, "neuron" = 10),
                                      algorithm.name = "nnls", return.dimensions = "tall")

dfp.wide <- get_ypb_experiment_series(sce.k2, sample.id.variable = "donor", 
                                      celltype.variable = "k2", assay.name = "logcounts",
                                      s.vector = c("glial" = 3, "neuron" = 10),
                                      algorithm.name = "nnls", return.dimensions = "wide")

# make new plots
# plot proportions panel -- no scale
ggplot(dfp.tall[dfp.tall$type=="noscale",], 
       aes(x = neuron.true, y = neuron.pred)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_hline(yintercept = 0.5) + geom_vline(xintercept = 0.5) + theme_bw() +
  xlab("True proportion") + ylab("Predicted proportion") +
  xlim(0, 1) + ylim(0, 1) + ggtitle("Neuron")

# plot proportions multipanel -- scale vs with scale
ggplot(dfp.tall, aes(x = neuron.true, y = neuron.pred)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_hline(yintercept = 0.5) + geom_vline(xintercept = 0.5) + theme_bw() +
  xlab("True proportion") + ylab("Predicted proportion") +
  xlim(0, 1) + ylim(0, 1) + facet_wrap(~type) + ggtitle("Neuron")

# jitterbox -- jittered points and boxplots of absolute errors
ggplot(dfp.tall, aes(x = type, y = neuron.abs.error)) + geom_jitter(alpha = 0.5) + 
  geom_boxplot(color = "cyan", alpha = 0) + theme_bw() + ggtitle("Neuron")

# save env
save.image("./cohort2/env/01_pseudobulk/01_k2.RData")
