#!/usr/bin/env R

#
# Test independent pseudobulk from DLPFC cohort1 snRNAseq data.
#

source("deconvo_method-paper/code/02_pseudobulk-predictions/00_parameters-pseudobulk.R")
sapply(libv, library, character.only = T)
list.sce.markers <- get(load(sce.markers.list.path))
sce.k2 <- list.sce.markers$k2
sce.k3 <- list.sce.markers$k3
sce.k4 <- list.sce.markers$k4

# get logcounts
sce.k2 <- logNormCounts(sce.k2)
sce.k3 <- logNormCounts(sce.k3)
sce.k4 <- logNormCounts(sce.k4)

# get pseudobulks
y.k2.unscale <- ypb_from_sce(sce = sce.k2, 
                             assay.name = "logcounts", 
                             celltype.variable = "k2", 
                             sample.id.variable = "Sample") %>% as.matrix()
y.k3.unscale <- ypb_from_sce(sce = sce.k3, 
                             assay.name = "logcounts", 
                             celltype.variable = "k3", 
                             sample.id.variable = "Sample") %>% as.matrix()
y.k4.unscale <- ypb_from_sce(sce = sce.k4, 
                             assay.name = "logcounts", 
                             celltype.variable = "k4", 
                             sample.id.variable = "Sample") %>% as.matrix()

# get true proportions
prop.true.k2 <- get_true_proportions(sce.k2, "k2", "donor")
prop.true.k3 <- get_true_proportions(sce.k3, "k3", "donor")
prop.true.k4 <- get_true_proportions(sce.k4, "k4", "donor")

# get predictions
filter.k2 <- rownames(sce.k2) %in% rownames(list.sce.markers$k2)
filter.k3 <- rownames(sce.k3) %in% rownames(list.sce.markers$k3)
filter.k4 <- rownames(sce.k4) %in% rownames(list.sce.markers$k4)
prop.pred.k2 <- lute(sce = sce.k2[filter.k2,], y = y.k2.unscale, 
                     assay.name = "logcounts", celltype.variable = "k2",
                     typemarker.algorithm = NULL)$deconvolution.results
prop.pred.k3 <- lute(sce = sce.k3[filter.k3,], y = y.k3.unscale,
                     assay.name = "logcounts", celltype.variable = "k3",
                     typemarker.algorithm = NULL)$deconvolution.results
prop.pred.k4 <- lute(sce = sce.k4[filter.k4,], y = y.k4.unscale,
                     assay.name = "logcounts", celltype.variable = "k4",
                     typemarker.algorithm = NULL)$deconvolution.results


