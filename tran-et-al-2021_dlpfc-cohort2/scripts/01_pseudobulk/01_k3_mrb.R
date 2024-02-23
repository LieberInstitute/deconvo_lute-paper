#!/usr/bin/env R

# Author: Sean Maden
#
# Test independent pseudobulk from DLPFC cohort2 snRNAseq data.
#
#
#
#

source("./cohort2/scripts/01_pseudobulk/00_parameters-pseudobulk.R")
sapply(libv, library, character.only = T)

# load marker data
markers.path <- "./cohort1/outputs/00_preprocess/list-scef_markers-k2-k3-k4_ro1-dlpfc.rda"
list.sce.markers <- get(load(markers.path))

# load sce data
sce.path <- "./cohort2/data/sce-mrb_dlpfc.rda"
sce.mrb <- get(load(sce.path)) # load sce data
sce.k3 <- sce.mrb[rownames(sce.mrb) %in% rownames(list.sce.markers[["k3"]]),]
s.vector <- c("Excit" = 10, "glial" = 3, "Inhib" = 10)


# get experiment results tables
dfp.tall <- get_ypb_experiment_series(sce.k3, sample.id.variable = "donor", 
                                      celltype.variable = "k3", assay.name = "logcounts",
                                      s.vector = s.vector,
                                      algorithm.name = "nnls", return.dimensions = "tall")
dfp.wide <- get_ypb_experiment_series(sce.k3, sample.id.variable = "donor", 
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

# jitterbox -- jittered points and boxplots of absolute errors
dfp.ae1 <- data.frame(celltype = dfp.ct$celltype,
                      abs.error = dfp.ct$abs.error.withscale)
dfp.ae1$type <- "withscale"
dfp.ae2 <- data.frame(celltype = dfp.ct$celltype,
                      abs.error = dfp.ct$abs.error.noscale)
dfp.ae2$type <- "noscale"
dfp.ae <- rbind(dfp.ae1, dfp.ae2)
ggplot(dfp.ae, aes(x = celltype, y = abs.error)) + geom_jitter(alpha = 0.5) + 
  geom_boxplot(color = "cyan", alpha = 0) + theme_bw() + facet_wrap(~type)



# save env
save.image("./cohort2/env/01_pseudobulk/01_k3_mrb_script.RData")
