#!/usr/bin/env R

# Author: Sean Maden
#
# Runs deconvolution with bias adjustment.
#

libv <- c("MultiAssayExperiment", "ggplot2", "gridExtra")
sapply(libv, library, character.only = T)

source("./source/00_dataset_summaries.R")
source("./source/00_deconvo_plots.R")
source("./source/00_sopt.R")
source("./source/00_musicParam-class.R")
source("./scripts/08_adjustment/00_param.R")

#-----
# load
#-----
mae <- get(load("./outputs/01_mae/mae_analysis_append.rda"))
dim(mae[[1]])

filter.mae <- colData(mae)$sample.id[complete.cases(mae)]
mae <- mae[,filter.mae,]
dim(mae[[1]])

sample.id.vector <- colData(mae)$sample.id
df.true.list <- metadata(mae[[1]])[["list.df.true.k2"]]
bulk.assay.name <- "bulk.pb.k2"
assay.name <- "counts"
y.unadj <- mae[[bulk.assay.name]]
dim(y.unadj)
sce <- mae[[1]]

#------------------------
# get s opt for 2 samples
#------------------------
# new run
dfs.param <- data.frame(
  s.min = c(1, NA, NA),
  s.max = c(80, NA, NA),
  s.step = c(5, 5e-2, 5e-3),
  s.diff = c(NA, 1, 1e-1)
)
list.res.all <- lapply(sample.id.vector, function(sample.id){
  list.res <- run_sopt_series(dfs.param, sample.id.vector = sample.id, 
                              df.true.list = df.true.list, y.unadj = y.unadj, 
                              sce = sce, assay.name = assay.name)
  list.res
})
names(list.res.all) <- sample.id.vector
min(list.res.all[[1]]$iter3$df.res$error.neuron)
min(list.res.all[[2]]$iter3$df.res$error.neuron)
list.dfsopt <- lapply(sample.id.vector, function(sample.id){
  dfs.opt <- list.res.all[[sample.id]][[nrow(dfs.param)]]$df.res
  dfs.opt <- dfs.opt[dfs.opt$error.glial==min(dfs.opt$error.glial),]
  matrix(c(sample.id, dfs.opt$s.glial[1], 
           dfs.opt$s.neuron[1], dfs.opt$error.neuron[1]), nrow = 1)
})
df.sopt <- do.call(rbind, lapply(list.dfsopt, function(item){item}))
df.sopt <- as.data.frame(df.sopt)
colnames(df.sopt) <- c("sample.id", "s.glial", "s.neuron", "min.error.neuron")
df.sopt

#-------------
# sanity check
#-------------
for(sample.id in sample.id.vector){
  mae.iter <- mae[,colData(mae)$sample.id==sample.id,]
  sce <- mae.iter[["snrnaseq.k2.all"]]
  sce <- scuttle::logNormCounts(sce)
  filt.dfs <- df.sopt$sample.id==sample.id
  dff.s <- df.sopt[filt.dfs,]
  s.vector.scale <- c("glial" = as.numeric(dff.s[2]), 
                      "neuron" = as.numeric(dff.s[3]))
  y.set <- mae.iter[[bulk.assay.name]]
  #y.set <- y.set[,y.set$expt_condition=="Nuc_RiboZeroGold"]
  y.set <- scuttle::logNormCounts(y.set)
  y <- assays(y.set)[[assay.name]][,,drop=F]
  nnls.noscale <- lute(sce, y = y, 
                       celltype.variable = "k2", 
                       s = s.vector.scale, assay.name = assay.name,
                       typemarker.algorithm = NULL)$deconvolution.results@predictions.table
  neuron.true <- as.numeric(metadata(sce)[["list.df.true.k2"]][[sample.id]]["neuron"])
  message(abs(neuron.true-nnls.noscale$neuron))
}

#--------------------------------
# get bias adj experiment results
#--------------------------------
list.s.vector.scale <- lapply(list.dfsopt, function(item){
  c("glial" = as.numeric(item[2]), 
    "neuron" = as.numeric(item[3]))
})
names(list.s.vector.scale) <- sample.id.vector

list.experiment.results <- experiment_all_samples(
  sample.id.vector, list.s.vector.scale = list.s.vector.scale,
  mae, bulk.mae.name = bulk.assay.name, assay.name = assay.name
)

df.res1 <- list.experiment.results$Br8667_ant$df.res

df.res <- do.call(rbind, lapply(list.experiment.results, function(item){
  df.res.iter <- item$df.res
  df.res.iter <- df.res.iter[df.res.iter$sample.id,]
  df.res.iter[1,,drop=F]
}))
df.res <- as.data.frame(df.res)

abs(df.res$neuron.true-df.res$neuron.music.scale)
abs(df.res$neuron.true-df.res$neuron.nnls.scale)
abs(df.res$neuron.true-df.res$neuron.bisque.scale)
abs(df.res$neuron.true-df.res$neuron.bisque.noscale)

#--------------
# get plot data
#--------------
dfp.wide <- df.res[,grepl("neuron", colnames(df.res))]

dfp.tall <- rbind(
  data.frame(value = df.res$neuron.nnls.scale, 
             variable = rep("nnls.scale", nrow(df.res)), 
             sample.id = df.res$sample.id),
  data.frame(value = df.res$neuron.nnls.noscale, 
             variable = rep("nnls.noscale", nrow(df.res)), 
             sample.id = df.res$sample.id),
  data.frame(value = df.res$neuron.music.scale, 
             variable = rep("music.scale", nrow(df.res)), 
             sample.id = df.res$sample.id),
  data.frame(value = df.res$neuron.music.noscale, 
             variable = rep("music.noscale", nrow(df.res)), 
             sample.id = df.res$sample.id),
  data.frame(value = df.res$neuron.bisque.scale, 
             variable = rep("bisque.scale", nrow(df.res)), 
             sample.id = df.res$sample.id),
  data.frame(value = df.res$neuron.bisque.noscale, 
             variable = rep("bisque.noscale", nrow(df.res)), 
             sample.id = df.res$sample.id)
)
dfp.tall <- as.data.frame(dfp.tall)
dfp.tall$true <- df.res$neuron.true
dfp.tall$error <- abs(dfp.tall$true-dfp.tall$value)
dfp.tall$cell.type <- "neuron"
dfp.tall$scale <- grepl(".*\\.scale$", dfp.tall$variable)
dfp.tall$algorithm <- gsub("\\..*", "", dfp.tall$variable)

list.dfp.plots <- get_dfp_list(df.res)


#-----
# save
#-----
# save s optima
save(df.sopt, file = "./outputs/09_fast/df_soptimize_pseudobulk_all.rda")

# save image
rm(mae)
save.image(file = "./env/09_fast/02_run_adjustment_pseudobulk_all_script.RData")
