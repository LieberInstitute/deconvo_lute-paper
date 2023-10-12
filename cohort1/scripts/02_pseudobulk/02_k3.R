#!/usr/bin/env R

# Author: Sean Maden
#
# Test independent pseudobulk from DLPFC cohort1 snRNAseq data.
#
#
#
#

source("scripts/02_pseudobulk/00_param.R")
sapply(libv, library, character.only = T)

# load mae samples
mae.filename <- "mae_analysis_append.rda"
mae.path <- file.path("outputs", "01_mae", mae.filename)
mae.analysis <- get(load(mae.path))
sce.k3 <- mae.analysis[["snrnaseq.k3.all"]]
s.vector <- c("Excit" = 10, "glial" = 3, "Inhib" = 10)

# get experiment results tables
dfp.tall <- get_ypb_experiment_series(sce.k3, sample.id.variable = "Sample", 
                                      celltype.variable = "k3", assay.name = "logcounts",
                                      s.vector = s.vector,
                                      algorithm.name = "nnls", return.dimensions = "tall")

dfp.wide <- get_ypb_experiment_series(sce.k3, sample.id.variable = "Sample", 
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



# save environment
save.image("./env/02_pseudobulk/02_k3.RData")