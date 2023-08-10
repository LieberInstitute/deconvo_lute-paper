#!/usr/bin/env R

# Author: Sean Maden
#
# Get full run of bias predictions.
#
# CRUCIAL NOTES, READ THIS:
#   * Z is the same across experiments
#   * Y is different across experiments
#

libv <- c("snow", "dplyr", "parallel", "doParallel", "lute", "dplyr")
sapply(libv, library, character.only = T)

#------
# setup
#------
# helper functions
# get series of s cell size factors (THIS SCRIPT, AND A FEW OTHERS)
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

# get bias computations in parallel (THIS SCRIPT, AND A FEW OTHERS)
parallel_bias_matched <- function(sce, yunadj, dfs, 
                                  celltype.variable = "k2",
                                  assay.name = "counts", 
                                  s.vector.ypb = c("glial" = 3, "neuron" = 10)){
  # begin parallel
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  # get full run
  if(is(yunadj, "RangedSummarizedExperiment")){
    yunadj <- assays(yunadj)[[1]] %>% as.matrix()
  }
  df.res <- do.call(rbind, 
                    mclapply(seq(nrow(dfs)), 
                             function(i){
                               s.vector <- c("glial" = dfs$glial[i], 
                                             "neuron" = dfs$neuron[i])
                               suppressMessages(
                                 lute(sce, y = yunadj, celltype.variable = celltype.variable, s = s.vector,
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

# get s vector series (THIS SCRIPT)
dfs <- dfs.series()

#----------
# load data
#----------
# load mae (SEE CODE 01 OUTPUTS)
new.mae.filename <- "mae_with-rpkm_additional-data_final.rda"
mae.final.filepath <- file.path("deconvo_method-paper", "outputs", "01_prepare-datasets", new.mae.filename)
mae <- get(load(mae.final.filepath))
sce <- mae[["sn1.rnaseq"]]
y.unadj <- mae[["bulk.rnaseq"]]

#------------
# main script
#------------
# set params (SEE PROJECT NOTES)
assay.name <- "counts"
celltype.variable <- "k2"
sample.id.variable <- "Sample"
sample.id.vector <- unique(sce[[sample.id.variable]])

# this is the chunk that makes the results df (CHECK CRUCIAL NOTES)
df.res.samples <- parallel_bias_matched(sce, y.unadj, dfs,
                                        celltype.variable = celltype.variable,
                                        assay.name = assay.name)


# append coldata from y.unadj (see MAE data)
df.res.samples$sample.labels <- rep(colnames(y.unadj), nrow(dfs))
df.res.samples$sample.id <- rep(gsub("_.*", "", y.unadj$batch.id), nrow(dfs))
df.res.samples$cell.compartment <- rep(y.unadj$expt_condition, nrow(dfs))
df.res.samples$anatomic.region <- rep(y.unadj$location, nrow(dfs))
df.res.samples$library.preparation <- rep(y.unadj$library_prep, nrow(dfs))
df.res.samples$sample.id.brnum <- rep(y.unadj$batch.id2, nrow(dfs))

# append data transformations
# this is the chunk that sets more operants in `df.res`
df.res.samples$s.fraction.neuron.glial <- df.res.samples$neuron/df.res.samples$glial
df.res.samples$log.s.fraction <- log(df.res.samples$s.fraction.neuron.glial)
df.res.samples$minimum.error <- df.res.samples$error.neuron==min(df.res.samples$minimum.error)
df.res.samples$maximum.error <- df.res.samples$error.neuron==max(df.res.samples$maximum.error)
deciles.error.neuron <- quantile(df.res.samples$error.neuron, seq(0, 1, 0.1))
df.res.samples$minimum.decile.error <- df.res.samples$error.neuron <= deciles.error.neuron[2]
df.res.samples$maximum.decile.error <- df.res.samples$error.neuron >= deciles.error.neuron[9]
df.res.samples$error.neuron <- df.res.samples$bias.neuron.true.pred %>% abs()

# save
save.filename <- "df-sopt-result_yvary-zsame_cohort1.rda"
save.path <- file.path("deconvo_method-paper", "outputs", 
                       "12_soptimize_yvary-zsame_dlpfc-cohort1", save.filename)
save(df.res.samples, file = save.path)
