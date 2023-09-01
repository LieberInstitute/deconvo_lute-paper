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
# get series of s cell size factors
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
# get bias computations in parallel
parallel_bias <- function(sce, dfs, celltype.variable = "k2", 
                          assay.name = "counts", s.vector.ypb = c("glial" = 3, "neuron" = 10)){
  # begin parallel
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  # get full run
  df.res <- do.call(rbind, 
                    mclapply(seq(nrow(dfs)), 
                             function(i){
                               s.vector <- c("glial" = dfs$glial[i], "neuron" = dfs$neuron[i])
                               ypb <- ypb_from_sce(sce, assay.name, 
                                                   celltype.variable, S = s.vector.ypb) %>% as.matrix()
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
  return(df.res)
}

# get s vector series
dfs <- dfs.series()

# load data
sce.markers.list.path <- file.path("deconvo_method-paper", "outputs", "01_prepare-datasets", 
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
sample.id.variable <- "Sample"
sample.id.vector <- unique(sce[[sample.id.variable]])

df.res.samples <- do.call(rbind, lapply(sample.id.vector, function(sample.id){
  message(sample.id)
  sce.iter <- sce[,sce[[sample.id.variable]]==sample.id]
  df.iter <- parallel_bias(sce.iter, dfs, 
                           celltype.variable = celltype.variable, 
                           assay.name = assay.name)
  df.iter$sample.id <- sample.id
  return(df.iter)
}))

# save
save.filename <- "df-result_s-opt-bias_cohort1.rda"
save.path <- file.path("deconvo_method-paper", "outputs", 
                       "11_soptimize-pbfit-bias_dlpfc-cohort1", save.filename)
save(df.res.samples, file = save.path)
