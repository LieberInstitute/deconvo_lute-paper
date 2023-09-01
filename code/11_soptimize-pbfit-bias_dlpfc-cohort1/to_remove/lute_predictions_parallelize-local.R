#!/usr/bin/env R

# Author: Sean Maden
#
# Time a limited run of nnls predictions using parallelization on Windows.
#

library(snow)

# helper functions
dfs.series <- function(
    s.glial.series = 
      seq(1, 20, 0.2)){
  s.neuron.series <- 
    rev(s.glial.series)
  dfs.series <- 
    do.call(rbind, 
            lapply(seq(length(s.glial.series)), function(index1){
    
              do.call(rbind, 
                      lapply(seq
                             (length
                               (s.neuron.series)), function(index2){
      
                                 c("glial" = s.glial.series[index1], 
                                   "neuron" = s.neuron.series[index2])
    
                                 }))
  
              })) %>% as.data.frame()
  
  return(dfs.series)
  
}

dfs <- dfs.series()

dim(dfs)

# load data
sce.markers.list.path <- file.path("deconvo_method-paper", 
                                   "outputs", "01_prepare-datasets", 
                                   "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda")
list.sce.markers <- get(load(sce.markers.list.path))
sce.k2 <- list.sce.markers$k2

# set params
sce <- sce.k2
assay.name <- "counts"
celltype.variable <- "k2"

# set dfs s variable
dfs <- dfs.series()

# time a limited run
cl <- makeCluster(detectCores())
registerDoParallel(cl)
system.time(
  mclapply(
    1:100, function(i){
    s.vector <- c("glial" = dfs$glial[i], "neuron" = dfs$neuron[i])
    ypb <- ypb_from_sce(
      sce, assay.name, celltype.variable) %>% as.matrix()
    suppressMessages(
      lute(
        sce.k2, y = ypb, celltype.variable = celltype.variable, 
                        typemarker.algorithm = NULL)$deconvolution.results@predictions.table)
}))
registerDoSEQ()