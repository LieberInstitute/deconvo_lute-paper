#!/usr/bin/env R

#
# Test independent pseudobulk from multi-region brain data.
#

source("deconvo_method-paper/code/02_pseudobulk-predictions/00_parameters-pseudobulk.R")
sapply(libv, library, character.only = T)
sce.mrb <- get(load(sce.mrb.path))
list.sce.markers <- get(load(sce.markers.list.path))

# filter on markers
all.markers <- unique(c(rownames(list.sce.markers[[1]]), 
                        rownames(list.sce.markers[[2]]),
                        rownames(list.sce.markers[[3]])))
sce.mrb <- sce.mrb[rownames(sce.mrb) %in% all.markers,]

# get pseudobulks
y.k2.unscale <- ypb_from_sce(sce = sce.mrb, 
                             assay.name = "logcounts", 
                             celltype.variable = "k2", 
                             sample.id.variable = "donor") %>% as.matrix()
y.k3.unscale <- ypb_from_sce(sce = sce.mrb, 
                             assay.name = "logcounts", 
                             celltype.variable = "k3", 
                             sample.id.variable = "donor") %>% as.matrix()
y.k4.unscale <- ypb_from_sce(sce = sce.mrb, 
                             assay.name = "logcounts", 
                             celltype.variable = "k4", 
                             sample.id.variable = "donor") %>% as.matrix()
 
# get true proportions
prop.true.k2 <- get_true_proportions(sce.mrb, "k2", "donor")
prop.true.k3 <- get_true_proportions(sce.mrb, "k3", "donor")
prop.true.k4 <- get_true_proportions(sce.mrb, "k4", "donor")

# get predictions
filter.k2 <- rownames(sce.mrb) %in% rownames(list.sce.markers$k2)
filter.k3 <- rownames(sce.mrb) %in% rownames(list.sce.markers$k3)
filter.k4 <- rownames(sce.mrb) %in% rownames(list.sce.markers$k4)
prop.pred.k2 <- lute(sce = sce.mrb[filter.k2,], y = y.k2.unscale, 
                        assay.name = "logcounts", celltype.variable = "k2",
                        typemarker.algorithm = NULL)$deconvolution.results
prop.pred.k3 <- lute(sce = sce.mrb[filter.k3,], y = y.k3.unscale,
                        assay.name = "logcounts", celltype.variable = "k3",
                        typemarker.algorithm = NULL)$deconvolution.results
prop.pred.k4 <- lute(sce = sce.mrb[filter.k4,], y = y.k4.unscale,
                        assay.name = "logcounts", celltype.variable = "k4",
                        typemarker.algorithm = NULL)$deconvolution.results

# get plot data
identical(rownames(prop.pred.k2), rownames(prop.true.k2))
identical(rownames(prop.pred.k3), rownames(prop.true.k3))
identical(rownames(prop.pred.k4), rownames(prop.true.k4))
colnames(prop.pred.k2) <- paste0(colnames(prop.pred.k2), ".pred")
colnames(prop.pred.k3) <- paste0(colnames(prop.pred.k3), ".pred")
colnames(prop.pred.k4) <- paste0(colnames(prop.pred.k4), ".pred")
colnames(prop.true.k2) <- paste0(colnames(prop.true.k2), ".true")
colnames(prop.true.k3) <- paste0(colnames(prop.true.k3), ".true")
colnames(prop.true.k4) <- paste0(colnames(prop.true.k4), ".true")
dfp.k2 <- cbind(prop.true.k2, prop.pred.k2) %>% as.data.frame()
dfp.k3 <- cbind(prop.true.k3, prop.pred.k3) %>% as.data.frame()
dfp.k4 <- cbind(prop.true.k4, prop.pred.k4) %>% as.data.frame()

# scatterplots of neurons
ggplot(dfp.k2, aes(x = neuron.true, y = neuron.pred)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0) + theme_bw() +
  xlab("True proportion") + ylab("Predicted proportion") +
  xlim(0, 1) + ylim(0, 1)
