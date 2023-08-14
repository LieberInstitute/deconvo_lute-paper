#!/usr/bin/env R

# Author: Sean Maden
# 
# Validate sopt utility on new held-out bulk samples.
#

libv <- c("lute", "MultiAssayExperiment")
sapply(libv, library, character.only = T)

# source
script.path <- file.path("deconvo_method-paper", "code", 
                         "12_soptimize_yvary-zsame_dlpfc-cohort1", "00_parameters.R")
source(script.path)

#-----
# load
#-----
# validation data
validation.sample.id <- c("Br6432", "Br6522", "Br8667")
# get all mae
mae.path <- file.path("deconvo_method-paper", "outputs", "01_prepare-datasets", "mae_with-rpkm_additional-data_final.rda")
mae <- get(load(mae.path))
# get sce marker data
sce <- mae[["sn1.rnaseq"]] # training marker expression
# subset mae
cd.mae <- colData(mae)
filter.validate.mae <- grepl(paste0(validation.sample.id,collapse = "|"),cd.mae$sample.id)
table(filter.validate.mae)
mae <- mae[,filter.validate.mae,]

# load the best performing s cell size scale factors
sopt.path <- file.path("deconvo_method-paper", "outputs", 
                       "12_soptimize_yvary-zsame_dlpfc-cohort1", 
                       "list-sopt-values_results-yvar_cohort1.rda")
list.sopt.train <- get(load(sopt.path))

#------------------------------------
# define the common z for experiments
#------------------------------------
sample.id <- "Br8492_post"
sce <- sce[,sce$Sample == sample.id]

#---------------------------------
# define the true cell proportions
#---------------------------------
df.rn <- mae[["df.cellstat.rnascope"]]
sample.id.vector <- validation.sample.id
list.df.true <- df.true.list(df.rn, sample.id.vector, "k2", c("glial", "neuron"))
names(list.df.true) <- sample.id.vector

#---------------------------
# get bulk validation subset
#---------------------------
y.validate <- mae[["bulk.rnaseq"]]
bulk.validate.filter <- y.validate$BrNum %in% validation.sample.id
y.validate <- y.validate[,bulk.validate.filter]
dim(y.validate)

#-----------
# assign dfs
#-----------
df.min.train <- list.sopt.train$minima.by.label
dfs <- data.frame(glial = df.min.train$glial, neuron = df.min.train$neuron)

#---------------------------
# get results for each s set
#---------------------------
# this is the chunk that makes the results df (CHECK CRUCIAL NOTES)
df.res.samples <- multigroup_bias_matched(validation.sample.id, 
                                          list.df.true, 
                                          y.validate, 
                                          dfs, 
                                          sce)


# set params (SEE PROJECT NOTES)
assay.name <- "counts"
celltype.variable <- "k2"
sample.id.variable <- "Sample"
sample.id.vector <- unique(sce[[sample.id.variable]])







