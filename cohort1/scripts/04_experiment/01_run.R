#!/usr/bin/env R

# Author: Sean Maden
#
# Run pseudobulk s-optimization experiments.
#

libv <- c("snow", "dplyr", "parallel", "doParallel", "lute", "dplyr")
sapply(libv, library, character.only = T)

# sets variables
folder.name <- "04_experiment"
assay.name <- "counts"

# source
script.path <- file.path("scripts", folder.name, "00_param.R")
source(script.path)

#----------
# load data
#----------
# load mae (SEE CODE 01 OUTPUTS)
new.mae.filename <- "mae_allsamples.rda"
mae.final.filepath <- file.path("outputs", "01_mae", new.mae.filename)
mae <- get(load(mae.final.filepath))

# remove validation samples
# dim(mae[["bulk.rnaseq"]])
cd.mae <- colData(mae)
cd.mae$sample.id.new <- gsub("_.*", "", cd.mae$sample.id)
validation.sample.id <- c("Br6432", "Br6522", "Br8667")
filter.string <- paste0(validation.sample.id, collapse = "|")
filter.mae <- !grepl(filter.string, cd.mae$sample.id.new)
table(filter.mae)
# filter.mae
# FALSE  TRUE 
# 7    15 

mae <- mae[,filter.mae,]
#dim(mae[["bulk.rnaseq"]])

# get assay data
sce <- mae[["snrnaseq.k2.all"]]
# y.unadj <- mae[["bulk.rnaseq"]]

#------------
# get y.unadj
#------------
sample.id.vector <- unique(sce[["Sample"]])
y.unadj.counts <- do.call(cbind, lapply(sample.id.vector, function(sample.id){
  ypb_from_sce(sce, "counts", "k2")
}))
colnames(y.unadj.counts) <- sample.id.vector

y.unadj <- SummarizedExperiment(assays = list(counts = as.matrix(y.unadj.counts)))
colData(y.unadj) <- DataFrame(data.frame(sample.id = sample.id.vector))
y.unadj <- as(y.unadj, "RangedSummarizedExperiment")
colnames(y.unadj) <- sample.id.vector

#----------------------------------
# get s vector series (THIS SCRIPT)
#----------------------------------
dfs <- dfs.series()

#---------------------------------
# define the true cell proportions
#---------------------------------
list.df.true <- lapply(sample.id.vector, function(sample.id){
  k.table <- table(sce[,sce[["Sample"]]==sample.id][["k2"]])
  k.prop <- as.data.frame(t(as.matrix(prop.table(k.table))))
  rownames(k.prop) <- "true_proportion"
  return(k.prop)
})
names(list.df.true) <- sample.id.vector

#------------
#
# main script
#
#------------
# set params (SEE PROJECT NOTES)
assay.name <- "counts"
celltype.variable <- "k2"
sample.id.variable <- "Sample"
sample.id.vector <- y.unadj$sample.id
df.res.samples <- do.call(rbind, lapply(sample.id.vector, function(sample.id){
  y.filt <- y.unadj[,y.unadj$sample.id==sample.id]
  sce.filt <- sce[,sce$Sample==sample.id]
  multigroup_bias_matched(sample.id, list.df.true, y.filt, 
                          y.group.name = "sample.id",
                          dfs, sce.filt, assay.name = assay.name)
}))
df.res.samples <- as.data.frame(df.res.samples)

# inspect
head(df.res.samples)
head(df.res.samples[df.res.samples$sample.label=="Br2720_mid",])

#--------------------
# postprocess results
#--------------------
# append data transformations
# this is the chunk that sets more operants in `df.res`
df.res.samples$s.fraction.neuron.glial <- df.res.samples$s.neuron/df.res.samples$s.glial
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
save.filename <- "pseudobulk-results.rda"
save.path <- file.path("deconvo_method-paper", "outputs", "04_experiment", save.filename)
save(df.res.samples, file = save.path)
