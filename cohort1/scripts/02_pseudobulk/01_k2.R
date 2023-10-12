#!/usr/bin/env R

# Author: Sean Maden
#
# Test independent pseudobulk from DLPFC cohort1 snRNAseq data.
#
#
#
#
#
#

source(file.path("scripts", "02_pseudobulk", "00_param.R"))
sapply(libv, library, character.only = T)

# load mae samples
mae.filename <- "mae_analysis_append.rda"
mae.path <- file.path("outputs", "01_mae", mae.filename)
mae.analysis <- get(load(mae.path))
sce.k2 <- mae.analysis[["snrnaseq.k2.all"]]

# get experiment results tables
dfp.tall <- get_ypb_experiment_series(sce.k2, sample.id.variable = "Sample", 
                                 celltype.variable = "k2", assay.name = "logcounts",
                                 s.vector = c("glial" = 3, "neuron" = 10),
                                 algorithm.name = "nnls", return.dimensions = "tall")

dfp.wide <- get_ypb_experiment_series(sce.k2, sample.id.variable = "Sample", 
                                      celltype.variable = "k2", assay.name = "logcounts",
                                      s.vector = c("glial" = 3, "neuron" = 10),
                                      algorithm.name = "nnls", return.dimensions = "wide")

dfp.ct <- dfp_tall_by_celltype(dfp.wide) # dfp.wide, tall by cell type


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

# save environment
save.image("./env/02_pseudobulk/01_k2.RData")

