#!/usr/bin/env R

# Author: Sean Maden
#
# Get full run of bias predictions.
#

library(snow)

#------
# setup
#------
# helper functions
dfs.series <- function(s.glial.series = seq(1, 20, 1)){
  s.neuron.series <- rev(s.glial.series)
  dfs.series <- do.call(rbind, lapply(seq(length(s.glial.series)), function(index1){
    do.call(rbind, lapply(seq(length(s.neuron.series)), function(index2){
      c("glial" = s.glial.series[index1], "neuron" = s.neuron.series[index2])
    }))
  })) %>% as.data.frame()
  #plot(dfs.series$glial, dfs.series$neuron)
  return(dfs.series)
}
dfs <- dfs.series()

# load data
sce.markers.list.path <- file.path("deconvo_method-paper", 
                                   "outputs", "01_prepare-datasets", 
                                   "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda")
list.sce.markers <- get(load(sce.markers.list.path))
sce.k2 <- list.sce.markers$k2

#------------
# main script
#------------
# set params
sce <- sce.k2
assay.name <- "counts"
celltype.variable <- "k2"

# begin parallel
cl <- makeCluster(detectCores())
registerDoParallel(cl)

# get full run
df.res <- do.call(rbind, mclapply(seq(nrow(dfs)), # seq(nrow(dfs)), 
                         function(i){
  s.vector <- c("glial" = dfs$glial[i], "neuron" = dfs$neuron[i])
  ypb <- ypb_from_sce(sce, assay.name, celltype.variable, S = s.vector) %>% as.matrix()
  suppressMessages(
    lute(sce, y = ypb, celltype.variable = celltype.variable, s = s.vector,
         typemarker.algorithm = NULL)$deconvolution.results@predictions.table
    )
}))
colnames(df.res) <- paste0(colnames(df.res), ".pred.nnls")
df.res <- cbind(df.res, dfs)
df.true <- sce[[celltype.variable]] %>% table() %>% prop.table() %>% as.data.frame()
rownames(df.true) <- df.true[,1]
df.res$glial.true <- df.true["glial",2]
df.res$neuron.true <- df.true["neuron",2]
df.res$bias.glial.true.pred <- df.res$glial.true - df.res$glial.pred.nnls
df.res$bias.neuron.true.pred <- df.res$neuron.true - df.res$neuron.pred.nnls

# make sequential again (i.e. cancels parallel)
registerDoSEQ()
