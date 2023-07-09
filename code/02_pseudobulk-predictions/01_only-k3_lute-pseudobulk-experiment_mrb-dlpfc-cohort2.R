#!/usr/bin/env R

#
# Test independent pseudobulk from DLPFC cohort1 snRNAseq data.
#

source("deconvo_method-paper/code/02_pseudobulk-predictions/00_parameters-pseudobulk.R")
sapply(libv, library, character.only = T)
list.sce.markers <- get(load(sce.markers.list.path))
sce <- list.sce.markers$k3
sce.mrb <- get(load(sce.mrb.path)) # load sce data
s.vector <- c("Excit" = 10, "glial" = 3, "Inhib" = 10)

# get experiment results tables
dfp.tall <- get_ypb_experiment_series(sce.mrb, sample.id.variable = "donor", 
                                      celltype.variable = "k3", assay.name = "logcounts",
                                      s.vector = s.vector,
                                      algorithm.name = "nnls", return.dimensions = "tall")
dfp.wide <- get_ypb_experiment_series(sce, sample.id.variable = "donor", 
                                      celltype.variable = "k3", assay.name = "logcounts",
                                      s.vector = s.vector,
                                      algorithm.name = "nnls", return.dimensions = "wide")
dfp.ct <- dfp_tall_by_celltype(dfp.wide) # dfp.wide, tall by cell type

# make new plots
# plot proportions multipanel
ggplot(dfp.ct, aes(x = true.noscale, y = pred.noscale)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_hline(yintercept = 0.5) + geom_vline(xintercept = 0.5) + theme_bw() +
  xlab("True proportion") + ylab("Predicted proportion") +
  xlim(0, 1) + ylim(0, 1) + facet_wrap(~celltype)
