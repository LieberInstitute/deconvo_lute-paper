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
# sample.id.remove <- c("Br2743_mid", "Br6432_mid", "Br6432_post")
filter.mae <- colData(mae)$sample.id[complete.cases(mae)]
mae <- mae[,filter.mae,]
dim(mae[[1]])

assay.name <- "logcounts"
bulk.assay.name <- "bulk.rnaseq"
y.unadj <- mae[[bulk.assay.name]]
y.unadj <- y.unadj[,y.unadj$expt_condition=="Nuc_RiboZeroGold"]
mae[[bulk.assay.name]] <- y.unadj
dim(y.unadj)
sce <- mae[[1]]
df.true.list <- metadata(sce)[["list.df.true.k2"]]
sample.id.vector <- unique(colData(mae)$sample.id[complete.cases(mae)])

#------------------------
# get s opt for 2 samples
#------------------------
# new run
dfs.param <- data.frame(
  s.min = c(1, NA, NA),
  s.max = c(80, NA, NA),
  s.step = c(2, 1e-2, 1e-3),
  s.diff = c(NA, 1, 1e-1)
)
list.res.all <- lapply(sample.id.vector, function(sample.id){
  y.unadj.iter <- y.unadj[,rep(which(y.unadj$batch.id2==sample.id),2)]
  sce.iter <- sce[,sce$Sample==sample.id]
  list.res <- run_sopt_series(dfs.param, sample.id.vector = sample.id, 
                              df.true.list = df.true.list[sample.id], 
                              y.unadj = y.unadj.iter, sce = sce.iter, 
                              assay.name = assay.name)
  list.res
})
names(list.res.all) <- sample.id.vector

for(sample.id in names(list.res.all)){
  message(sample.id)
  message(identical(list.res.all[[sample.id]][[3]]$df.res$neuron.true[1],
                    as.numeric(df.true.list[[sample.id]]["neuron"])))
}

min(list.res.all[[1]]$iter3$df.res$error.neuron)
min(list.res.all[[2]]$iter3$df.res$error.neuron)
list.dfsopt <- lapply(sample.id.vector, function(sample.id){
  dfs.opt <- list.res.all[[sample.id]][[nrow(dfs.param)]]$df.res
  dfs.opt <- dfs.opt[dfs.opt$error.neuron==min(dfs.opt$error.neuron),]
  matrix(c(sample.id, dfs.opt$s.glial[1], dfs.opt$s.neuron[1], 
           dfs.opt$error.neuron[1]), nrow = 1)
})
names(list.dfsopt) <- sample.id.vector
df.sopt <- do.call(rbind, lapply(list.dfsopt, function(item){item}))
df.sopt <- as.data.frame(df.sopt)
colnames(df.sopt) <- c("sample.id", "s.glial", "s.neuron", "min.error.neuron")
df.sopt

# exclude failed runs
max.allowed.error <- 0.1
df.sopt$min.error.neuron <- as.numeric(df.sopt$min.error.neuron)
filter.sopt <- df.sopt$min.error.neuron > 0.1
df.sopt[filter.sopt,]$sample.id # "Br6423_ant"  "Br6522_post"
df.sopt <- df.sopt[!filter.sopt,]
# update sample.id.vector
sample.id.vector <- df.sopt$sample.id

#-------------
# sanity check
#-------------
df.sopt$sanity.check.error <- NA
for(sample.id in sample.id.vector){
  mae.iter <- mae[,colData(mae)$sample.id==sample.id,]
  sce <- mae.iter[["snrnaseq.k2.all"]]
  sce <- scuttle::logNormCounts(sce)
  filt.dfs <- df.sopt$sample.id==sample.id
  dff.s <- df.sopt[filt.dfs,]
  s.vector.scale <- c("glial" = as.numeric(dff.s[2]), 
                      "neuron" = as.numeric(dff.s[3]))
  y.set <- mae.iter[[bulk.assay.name]]
  y.set <- y.set[,y.set$expt_condition=="Nuc_RiboZeroGold"]
  y.set <- scuttle::logNormCounts(y.set)
  y <- assays(y.set)[["logcounts"]][,,drop=F]
  nnls.noscale <- lute(sce, y = y, celltype.variable = "k2", 
                       s = s.vector.scale, assay.name = assay.name,
                       typemarker.algorithm = NULL)$deconvolution.results@predictions.table
  neuron.true <- as.numeric(metadata(sce)[["list.df.true.k2"]][[sample.id]]["neuron"])
  error.new <- abs(neuron.true-nnls.noscale$neuron)
  df.sopt[df.sopt$sample.id==sample.id,]$sanity.check.error <- error.new
}

df.sopt$min.error.neuron <- as.numeric(df.sopt$min.error.neuron)

ggplot(df.sopt, aes(x = min.error.neuron, y = sanity.check.error)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0)


#--------------------------------
# get bias adj experiment results
#--------------------------------
# filter samples that fail for bisque
#
# e.g. throws error:
# > Error in solve.default(t(D.weight) %*% D.weight) : 
# > Lapack routine dgesv: system is exactly singular: U[1,1] = 0
# 
filter.bisque.samples <- c("Br8667_ant")
sample.id.vector <- sample.id.vector[!sample.id.vector %in% filter.bisque.samples]

list.dfsopt <- list.dfsopt[names(list.dfsopt) %in% sample.id.vector]

list.s.vector.scale <- lapply(list.dfsopt, function(item){
  c("glial" = as.numeric(item[2]), "neuron" = as.numeric(item[3]))
})
names(list.s.vector.scale) <- sample.id.vector

list.experiment.results <- experiment_all_samples(
  sample.id.vector, list.s.vector.scale = list.s.vector.scale,
  mae, bulk.mae.name = "bulk.pb.k2", assay.name = "logcounts"
)

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
library(GGally)

dfp.wide <- df.res[,grepl("neuron", colnames(df.res))]
ggpairs(dfp.wide)

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

# replace missmatched true neuron proportions
for(sample.id in unique(dfp.tall$sample.id)){
  filter.dfp.tall <- dfp.tall$sample.id==sample.id
  filter.dfp.wide <- rownames(dfp.wide)==sample.id
  true.neuron.value <- as.numeric(df.true.list[[sample.id]]["neuron"])
  dfp.tall[filter.dfp.tall,]$true <- true.neuron.value
  dfp.wide[filter.dfp.wide,]$neuron.true <- true.neuron.value
}
dfp.tall$error <- abs(dfp.tall$true-dfp.tall$value)

#-----
# save
#-----
# save s optima
save(df.sopt, file = "./outputs/09_fast/df_soptimize_realbulk_all.rda")

# save image
rm(mae)
save.image(file = "./env/09_fast/03_run_adjustment_realbulk_all_script.RData")
