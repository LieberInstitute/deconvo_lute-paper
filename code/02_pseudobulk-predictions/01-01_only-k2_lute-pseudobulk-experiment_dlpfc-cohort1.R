#!/usr/bin/env R

#
# Test independent pseudobulk from DLPFC cohort1 snRNAseq data.
#

source("deconvo_method-paper/code/02_pseudobulk-predictions/00_parameters-pseudobulk.R")
sapply(libv, library, character.only = T)
list.sce.markers <- get(load(sce.markers.list.path))
sce.k2 <- list.sce.markers$k2

# get experiment results tables
dfp.tall <- get_ypb_experiment_series(sce.k2, sample.id.variable = "Sample", 
                                 celltype.variable = "k2", assay.name = "logcounts",
                                 s.vector = c("glial" = 3, "neuron" = 10),
                                 algorithm.name = "nnls", return.dimensions = "tall")

dfp.wide <- get_ypb_experiment_series(sce.k2, sample.id.variable = "Sample", 
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

# plot jitterbox absolute errors -- scale vs with scale
dfp.ae <- data.frame(abs.error.neuron.noscale = abs(dfp.wide$neuron.pred.noscale-dfp.wide$neuron.true.withscale),
                     abs.error.neuron.scale = abs(dfp.wide$neuron.pred.withscale-dfp.wide$neuron.true.withscale))
ggplot(dfp.tall, aes(x = type, y = abs.error.neuron)) + geom_jitter(alpha = 0.5) + 
  geom_boxplot(draw_quantiles = 0.5, alpha = 0, color = "cyan") + theme_bw() + facet_wrap(~type)

