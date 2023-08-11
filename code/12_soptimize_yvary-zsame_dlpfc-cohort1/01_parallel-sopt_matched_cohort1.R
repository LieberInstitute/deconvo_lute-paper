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

# parallel_bias_matched
# get bias computations in parallel (THIS SCRIPT, AND A FEW OTHERS)
parallel_bias_matched <- function(sce, yunadj, dfs, df.true = NULL, 
                                  celltype.variable = "k2", assay.name = "counts", 
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
                               s.vector <- c("glial" = dfs$glial[i], "neuron" = dfs$neuron[i])
                               suppressMessages(
                                 lute(sce, y = yunadj, celltype.variable = celltype.variable, s = s.vector,
                                      typemarker.algorithm = NULL)$deconvolution.results@predictions.table
                               )
                             }))
  colnames(df.res) <- paste0(colnames(df.res), ".pred.nnls")
  df.res <- cbind(df.res, dfs)
  if(is(df.true, "NULL")){
    df.true <- sce[[celltype.variable]] %>% table() %>% prop.table() %>% as.data.frame()
    rownames(df.true) <- df.true[,1]
  }
  df.res$glial.true <- as.numeric(df.true["glial",1])
  df.res$neuron.true <- as.numeric(df.true["neuron",1])
  df.res$bias.glial.true.pred <- df.res$glial.true - df.res$glial.pred.nnls
  df.res$bias.neuron.true.pred <- df.res$neuron.true - df.res$neuron.pred.nnls
  # make sequential again (i.e. cancels parallel)
  registerDoSEQ()
  return(df.res)
}

# multigroup_bias_matched
# wraps parallel_bias_matched for multiple groups, uses df.true.list
multigroup_bias_matched <- function(sample.id.vector, df.true.list, y.unadj, dfs, sce, 
                                    y.group.name = "batch.id2",
                                    celltype.variable = "k2", assay.name = "counts", 
                                    s.vector.ypb = c("glial" = 3, "neuron" = 10)){
  df.res <- do.call(rbind, lapply(sample.id.vector, function(sample.id){
    y.iter <- y.unadj[,colData(y.unadj)[,y.group.name]==sample.id]
    df.true.iter <- df.true.list[[sample.id]] %>% t() %>% as.data.frame()
    #colnames(df.true.iter.transpose) <- colnames(df.true.list[[sample.id]])
    df.res.iter <- parallel_bias_matched(sce, y.iter, dfs, df.true.iter, 
                                         celltype.variable, assay.name, s.vector.ypb)
    df.res.iter$sample.id <- sample.id
    return(df.res.iter)
  }))
  return(df.res)
}


# get true cell proportions info
rnascope_cell_info <- function(df.rn, sample.id = NULL, k.type = NULL, cell.types = NULL){
  #
  #
  # get cell info data from summary df
  #
  # example:
  # df.rn <- mae[["df.cellstat.rnascope"]]
  # ct.info <- rnascope_cell_info(df.rn, "k2", c("glial", "neuron"))
  #
  filter.vector <- colnames(df.rn)
  filter.condition.vector <- rep(TRUE, length(filter.vector))
  ct.string <- paste0(cell.types, collapse = "|")
  if(!is(ct.string, "NULL")){
    filter.condition.vector <- filter.condition.vector & grepl(ct.string, filter.vector)
  }
  if(!is(sample.id, "NULL")){
    filter.condition.vector <- filter.condition.vector & grepl(sample.id, filter.vector)
  }
  if(!is(k.type, "NULL")){
    filter.condition.vector <- filter.condition.vector & grepl(k.type, filter.vector)
  }
  df.ct <- df.rn[,filter.condition.vector]
  return(df.ct)
}
df.true.list <- function(df.rn, sample.id.vector, k.type, cell.types, info = "true_proportion"){
  #
  #
  # example:
  # sample.id.vector <- unique(y.unadj$batch.id2)
  # list.dftrue <- df.true.list(df.rn, sample.id.vector, "k2", c("glial", "neuron"))
  # names(list.dftrue) <- sample.id.vector
  #
  # NOTE: does not check cell type label order .. expect 1. glial, 2. neuron !!!
  #
  
  list.dfinfo <- lapply(sample.id.vector, function(sample.id){
    rnascope_cell_info(df.rn, sample.id, k.type, cell.types)
  })
  list.dfinfo <- lapply(list.dfinfo, function(dfinfo){
    dfinfo <- as.data.frame(dfinfo["true_proportion",,drop=F])
    colnames(dfinfo) <- cell.types
    return(dfinfo)
  })
  names(list.dfinfo) <- sample.id.vector
  return(list.dfinfo)
}

#----------
# load data
#----------
# load mae (SEE CODE 01 OUTPUTS)
new.mae.filename <- "mae_with-rpkm_additional-data_final.rda"
mae.final.filepath <- file.path("deconvo_method-paper", "outputs", "01_prepare-datasets", new.mae.filename)
mae <- get(load(mae.final.filepath))
sce <- mae[["sn1.rnaseq"]]
y.unadj <- mae[["bulk.rnaseq"]]

#----------------------------------
# get s vector series (THIS SCRIPT)
#----------------------------------
dfs <- dfs.series()

#---------------------------------
# define the true cell proportions
#---------------------------------
df.rn <- mae[["df.cellstat.rnascope"]]
sample.id.vector <- unique(y.unadj$batch.id2)
list.df.true <- df.true.list(df.rn, sample.id.vector, "k2", c("glial", "neuron"))
names(list.df.true) <- sample.id.vector

#------------------------------------
# define the common z for experiments
#------------------------------------
sample.id <- "Br8492_post"
sce <- sce[,sce$Sample == sample.id]

#------------
# main script
#------------
# set params (SEE PROJECT NOTES)
assay.name <- "counts"
celltype.variable <- "k2"
sample.id.variable <- "Sample"
sample.id.vector <- unique(sce[[sample.id.variable]])

# this is the chunk that makes the results df (CHECK CRUCIAL NOTES)
sample.id.vector <- unique(y.unadj$batch.id2)
df.res.samples <- multigroup_bias_matched(sample.id.vector, list.df.true, y.unadj, dfs, sce)

#df.res.samples <- parallel_bias_matched(sce, y.unadj, dfs, df.true,
#                                        celltype.variable = celltype.variable,
#                                        assay.name = assay.name)


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
df.res.samples$error.neuron <- abs(df.res.samples$bias.neuron.true.pred)
df.res.samples$error.glial <- abs(df.res.samples$bias.glial.true.pred)
df.res.samples$minimum.error <- df.res.samples$error.neuron==min(df.res.samples$error.neuron)
df.res.samples$maximum.error <- df.res.samples$error.neuron==max(df.res.samples$error.neuron)
deciles.error.neuron <- quantile(df.res.samples$error.neuron, seq(0, 1, 0.1))
df.res.samples$minimum.decile.error <- df.res.samples$error.neuron <= deciles.error.neuron[2]
df.res.samples$maximum.decile.error <- df.res.samples$error.neuron >= deciles.error.neuron[9]
df.res.samples$error.neuron <- df.res.samples$bias.neuron.true.pred %>% abs()

# save
save.filename <- "df-sopt-result_yvary-zsame_cohort1.rda"
save.path <- file.path("deconvo_method-paper", "outputs", 
                       "12_soptimize_yvary-zsame_dlpfc-cohort1", save.filename)
save(df.res.samples, file = save.path)
