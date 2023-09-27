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
sample.id.vector <- c("Br8667_ant", "Br8667_mid")
mae <- mae[,colData(mae)$sample.id %in% sample.id.vector,]
df.true.list <- metadata(mae[[1]])[["list.df.true.k2"]]
y.unadj <- mae[["bulk.rnaseq"]]
y.unadj <- y.unadj[,y.unadj$expt_condition=="Nuc_RiboZeroGold"]
mae[["bulk.rnaseq"]] <- y.unadj
dim(y.unadj)
sce <- mae[[1]]

#------------------------
# get s opt for 2 samples
#------------------------
# new run
dfs.param <- data.frame(
  s.min = c(1, NA, NA),
  s.max = c(40, NA, NA),
  s.step = c(1, 1e-2, 1e-3),
  s.diff = c(NA, 1, 1e-1)
)
list.res.all <- lapply(sample.id.vector, function(sample.id){
  sce <- scuttle::logNormCounts(sce)
  y.unadj <- scuttle::logNormCounts(y.unadj)
  list.res <- run_sopt_series(dfs.param, 
                              sample.id.vector = sample.id, 
                              df.true.list = df.true.list, 
                              y.unadj = y.unadj, 
                              sce = sce, 
                              assay.name = "logcounts")
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
# sample.id s.glial s.neuron     min.error.neuron
# 1 Br8667_ant  26.981    1.981 3.25372717460692e-05
# 2 Br8667_mid  37.013    3.013 1.05619995389716e-05

#-------------
# sanity check
#-------------
sample.id <- "Br8667_ant"
mae.iter <- mae[,colData(mae)$sample.id==sample.id,]
sce <- mae.iter[["snrnaseq.k2.all"]]
s.vector.scale <- c("glial" = as.numeric(df.sopt[1,2]), 
                    "neuron" = as.numeric(df.sopt[1,3]))
y.set <- mae.iter[["bulk.rnaseq"]]
y <- assays(y.set)[["counts"]][,,drop=F]
sce <- scuttle::logNormCounts(sce)
z <- lute::get_z_from_sce(sce, "counts", "k2")
nnls.noscale <- lute(z = z, y = y,
                     s = s.vector.scale, 
                     assay.name = "counts",
                     typemarker.algorithm = NULL)$deconvolution.results@predictions.table
nnls.noscale
#glial    neuron
#1 0.4994958 0.5005042
neuron.true <- as.numeric(metadata(sce)[["list.df.true.k2"]][[sample.id]]["neuron"])
abs(neuron.true-nnls.noscale$neuron)

sample.id <- "Br8667_mid"
mae.iter <- mae[,colData(mae)$sample.id==sample.id,]
sce <- mae.iter[["snrnaseq.k2.all"]]
s.vector.scale <- c("glial" = as.numeric(df.sopt[1,2]), 
                    "neuron" = as.numeric(df.sopt[1,3]))
y.set <- mae.iter[["bulk.rnaseq"]]
y <- assays(y.set)[["counts"]][,,drop=F]
sce <- scuttle::logNormCounts(sce)
z <- lute::get_z_from_sce(sce, "counts", "k2")
nnls.noscale <- lute(z = z, y = y,
                     s = s.vector.scale, 
                     assay.name = "counts",
                     typemarker.algorithm = NULL)$deconvolution.results@predictions.table
nnls.noscale
#glial    neuron
#1 0.4994958 0.5005042
neuron.true <- as.numeric(metadata(sce)[["list.df.true.k2"]][[sample.id]]["neuron"])
abs(neuron.true-nnls.noscale$neuron)

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
  mae, bulk.mae.name = "bulk.pb.k2"
)

df.res1 <- list.experiment.results$Br8667_ant$df.res
abs(df.res1$neuron.true-df.res1$neuron.music.scale)
abs(df.res1$neuron.true-df.res1$neuron.nnls.scale)

#-------
# save
#-------
# save image
rm(mae)
save.image(file = "./env/09_fast/02_run_adjustment_pseudobulk_2samples_script.RData")
