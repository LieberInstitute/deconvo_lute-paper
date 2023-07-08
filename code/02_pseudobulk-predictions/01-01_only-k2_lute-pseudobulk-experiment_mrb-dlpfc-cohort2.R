#!/usr/bin/env R

#
# Test independent pseudobulk from DLPFC cohort1 snRNAseq data.
#

source("deconvo_method-paper/code/02_pseudobulk-predictions/00_parameters-pseudobulk.R")
sapply(libv, library, character.only = T)
list.sce.markers <- get(load(sce.markers.list.path))

list.sce.markers <- get(load(sce.markers.list.path)) # load marker data
sce.mrb <- get(load(sce.mrb.path)) # load sce data
sce.k2 <- sce.mrb[rownames(sce.mrb) %in% rownames(list.sce.markers[["k2"]]),]
dfp.k2 <- get_ypb_experiment_results(sce.k2, sample.id.variable = "donor", 
                                     celltype.variable = "k2", assay.name = "logcounts",
                                     s.vector = c("glial" = 3, "neuron" = 10))

# scatterplots of neurons
ggplot(dfp.k2, aes(x = neuron.true, y = neuron.pred)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_hline(yintercept = 0.5) + geom_vline(xintercept = 0.5) + theme_bw() +
  xlab("True proportion") + ylab("Predicted proportion") +
  xlim(0, 1) + ylim(0, 1)
