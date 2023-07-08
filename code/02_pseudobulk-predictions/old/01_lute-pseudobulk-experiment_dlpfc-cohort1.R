#!/usr/bin/env R

#
# Test independent pseudobulk from DLPFC cohort1 snRNAseq data.
#

source("deconvo_method-paper/code/02_pseudobulk-predictions/00_parameters-pseudobulk.R")
sapply(libv, library, character.only = T)
list.sce.markers <- get(load(sce.markers.list.path))
sce.k2 <- list.sce.markers$k2
sce.k3 <- list.sce.markers$k3

# get logcounts
sce.k2 <- logNormCounts(sce.k2)
sce.k3 <- logNormCounts(sce.k3)

# get pseudobulks
y.k2.unscale <- ypb_from_sce(sce = sce.k2, 
                             assay.name = "logcounts", 
                             celltype.variable = "k2", 
                             sample.id.variable = "Sample") %>% as.matrix()
y.k3.unscale <- ypb_from_sce(sce = sce.k3, 
                             assay.name = "logcounts", 
                             celltype.variable = "k3", 
                             sample.id.variable = "Sample") %>% as.matrix()

# get true proportions
prop.true.k2 <- get_true_proportions(sce.k2, "k2", "Sample")
prop.true.k3 <- get_true_proportions(sce.k3, "k3", "Sample")

# get predictions
filter.k2 <- rownames(sce.k2) %in% rownames(list.sce.markers$k2)
filter.k3 <- rownames(sce.k3) %in% rownames(list.sce.markers$k3)
prop.pred.k2 <- lute(sce = sce.k2[filter.k2,], y = y.k2.unscale, 
                     assay.name = "logcounts", celltype.variable = "k2",
                     typemarker.algorithm = NULL)$deconvolution.results
prop.pred.k3 <- lute(sce = sce.k3[filter.k3,], y = y.k3.unscale,
                     assay.name = "logcounts", celltype.variable = "k3",
                     typemarker.algorithm = NULL)$deconvolution.results

# get plot data
identical(rownames(prop.pred.k2), rownames(prop.true.k2))
identical(rownames(prop.pred.k3), rownames(prop.true.k3))
colnames(prop.pred.k2) <- paste0(colnames(prop.pred.k2), ".pred")
colnames(prop.pred.k3) <- paste0(colnames(prop.pred.k3), ".pred")
colnames(prop.true.k2) <- paste0(colnames(prop.true.k2), ".true")
colnames(prop.true.k3) <- paste0(colnames(prop.true.k3), ".true")
dfp.k2 <- cbind(prop.true.k2, prop.pred.k2) %>% as.data.frame()
dfp.k3 <- cbind(prop.true.k3, prop.pred.k3) %>% as.data.frame()

# scatterplots of neurons
ggplot(dfp.k2, aes(x = neuron.true, y = neuron.pred)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0) + theme_bw() +
  xlab("True proportion") + ylab("Predicted proportion") +
  xlim(0, 1) + ylim(0, 1)

