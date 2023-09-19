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

# sets variables
folder.name <- "04_experiment_bulk"
assay.name <- "logcounts"

# source
script.path <- file.path("scripts", folder.name, "00_parameters.R")
source(script.path)

#----------
# load data
#----------
# load mae (SEE CODE 01 OUTPUTS)
new.mae.filename <- "mae_allsamples.rda"
mae.final.filepath <- file.path("outputs", "01_mae", new.mae.filename)
mae <- get(load(mae.final.filepath))

# remove validation samples
dim(mae[["bulk.rnaseq"]])
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
dim(mae[["bulk.rnaseq"]])

# get assay data
sce <- mae[["snrnaseq.k2.all"]]
y.unadj <- mae[["bulk.rnaseq"]]

#----------------------------------
# get s vector series (THIS SCRIPT)
#----------------------------------
dfs <- dfs.series()

#---------------------------------
# define the true cell proportions
#---------------------------------
df.rn <- mae[["cell.sizes"]]
sample.id.vector <- unique(y.unadj$batch.id2)
list.df.true <- df.true.list(df.rn, sample.id.vector, "k2", c("glial", "neuron"))
names(list.df.true) <- sample.id.vector

#------------------------------------
# define the common z for experiments
#------------------------------------
sample.id <- "Br8492_post"
sce <- sce[,sce$Sample == sample.id]

#------------
#
# main script
#
#------------
# set params (SEE PROJECT NOTES)
assay.name <- "counts"
celltype.variable <- "k2"
sample.id.variable <- "Sample"
sample.id.vector <- unique(sce[[sample.id.variable]])

# this is the chunk that makes the results df (CHECK CRUCIAL NOTES)
sample.id.vector <- unique(y.unadj$batch.id2)
df.res.samples <- multigroup_bias_matched(sample.id.vector, list.df.true, y.unadj, dfs, sce, assay.name = assay.name)
# inspect
head(df.res.samples)
head(df.res.samples[df.res.samples$sample.label=="2107UNHS-0291_Br2720_Mid_Bulk",])

#--------------------
# postprocess results
#--------------------
# append coldata from y.unadj (see MAE data)
cd.ydata <- colData(y.unadj)
#df.res.samples$sample.labels <- cd.ydata[,] rep(colnames(y.unadj), nrow(dfs))
#df.res.samples$sample.id <- rep(gsub("_.*", "", y.unadj$batch.id), nrow(dfs))
df.res.samples$cell.compartment <- cd.ydata[df.res.samples$sample.label,]$library_prep
df.res.samples$anatomic.region <- cd.ydata[df.res.samples$sample.label,]$location
df.res.samples$library.type <- cd.ydata[df.res.samples$sample.label,]$library_type
df.res.samples$sample.id.brnum <- cd.ydata[df.res.samples$sample.label,]$batch.id2

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
save.filename <- "df-sopt-result_yvary-zsame_cohort1.rda"
save.path <- file.path("deconvo_method-paper", "outputs", folder.name, save.filename)
save(df.res.samples, file = save.path)
