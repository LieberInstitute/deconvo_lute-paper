#!/usr/bin/env R

#
# Test independent pseudobulk from DLPFC cohort1 snRNAseq data.
#

source("deconvo_method-paper/code/02_pseudobulk-predictions/00_parameters-pseudobulk.R")
sapply(libv, library, character.only = T)
list.sce.markers <- get(load(sce.markers.list.path))
sce.k2 <- list.sce.markers$k2

get_ypb_experiment_results <- function(sce, sample.id.variable = "Sample", 
                                       celltype.variable = "k2", assay.name = "logcounts",
                                       s.vector = c("glial" = 3, "neuron" = 10)){
  unique.sample.id.vector <- unique(sce[[sample.id.variable]])
  dfp <- do.call(rbind, lapply(unique.sample.id.vector, function(sample.id){
    sce.iter <- sce[,sce$Sample==sample.id]
    ypb.iter <- ypb_from_sce(sce = sce.iter, assay.name = assay.name, 
                             celltype.variable = celltype.variable,
                             sample.id.variable = sample.id.variable, 
                             S = s.vector) %>% as.matrix()
    prop.true.iter <- table(sce.iter[[celltype.variable]]) %>% prop.table() %>% as.matrix() %>% t()
    prop.pred.iter <- lute(sce = sce.iter, y = ypb.iter, assay.name = assay.name, 
                           celltype.variable = celltype.variable, typemarker.algorithm = NULL, 
                           return.info = FALSE)$deconvolution.results
    colnames(prop.pred.iter) <- paste0(colnames(prop.pred.iter), ".pred")
    colnames(prop.true.iter) <- paste0(colnames(prop.true.iter), ".true")
    dfp.iter <- cbind(prop.true.iter, prop.pred.iter) %>% as.data.frame()
    dfp.iter
  }))
  rownames(dfp) <- unique.sample.id.vector
  return(dfp)
}

dfp.k2 <- get_ypb_experiment_results(sce.k2)

# scatterplots of neurons
ggplot(dfp.k2, aes(x = neuron.true, y = neuron.pred)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_hline(yintercept = 0.5) + geom_vline(xintercept = 0.5) + theme_bw() +
  xlab("True proportion") + ylab("Predicted proportion") +
  xlim(0, 1) + ylim(0, 1)
