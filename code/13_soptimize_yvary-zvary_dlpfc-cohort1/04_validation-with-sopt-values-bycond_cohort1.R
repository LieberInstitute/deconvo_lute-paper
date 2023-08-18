#!/usr/bin/env R

# Author: Sean Maden
# 
# Validate sopt utility on new held-out bulk samples.
#

libv <- c("lute", "MultiAssayExperiment")
sapply(libv, library, character.only = T)

# source
folder.name <- "13_soptimize_yvary-zvary_dlpfc-cohort1"
script.path <- file.path("deconvo_method-paper", "code", folder.name, "00_parameters.R")
source(script.path)

#-----
# load
#-----
# get all mae
mae.path <- file.path("deconvo_method-paper", "outputs", "01_prepare-datasets", "mae_final.rda")
mae <- get(load(mae.path))

# get sce marker data
# sce <- mae[["sn1.validate.rnaseq"]] # training marker expression
sce <- get(
  load(
    file.path("deconvo_method-paper", "outputs", "01_prepare-datasets", "list-sce-validate_markers-k-2-3-4_cohort1.rda")))[["k2"]]
sce <- sce[,!colData(sce)$k2=="other"]

# get df.rn from mae
mae.dfrn.path <- file.path("deconvo_method-paper", "outputs", "01_prepare-datasets", "mae_with-rpkm_additional-data_final.rda")
mae.dfrn <- get(load(mae.dfrn.path))
df.rn <- mae.dfrn[["df.cellstat.rnascope"]]
rm(mae.dfrn)

#-----------------------
# filter validation data
#-----------------------
# validation data
validation.sample.id.vector <- c("Br6432", "Br6522", "Br8667")
# subset mae
cd.mae <- colData(mae)
filter.validate.mae <- grepl(paste0(validation.sample.id.vector,collapse = "|"),cd.mae$sample.id)
table(filter.validate.mae)
mae <- mae[,filter.validate.mae,]
# load the best performing s cell size scale factors
sopt.path <- file.path("deconvo_method-paper", "outputs", "12_soptimize_yvary-zsame_dlpfc-cohort1", 
                       "list-sopt-values_results-yvar_cohort1.rda")
list.sopt.train <- get(load(sopt.path))

#---------------------------
# get bulk validation subset
#---------------------------
y.validate <- mae[["bulk.rnaseq"]]
rownames(y.validate) <- rowData(y.validate)$Symbol
bulk.validate.filter <- y.validate$BrNum %in% validation.sample.id.vector
y.validate <- y.validate[,bulk.validate.filter]
y.validate <- y.validate[rownames(sce),]
dim(y.validate)
# [1] 40 36

#---------------------------------
# define the true cell proportions
#---------------------------------
# df.rn <- mae[["df.cellstat.rnascope"]]
validation.sample.id.vector <- unique(y.validate$batch.id2)
list.df.true <- df.true.list(df.rn, validation.sample.id.vector, "k2", c("glial", "neuron"))
names(list.df.true) <- validation.sample.id.vector

#-----------
# assign dfs
#-----------
# load 
dfs.name <- "dfs-medians-bygroup-training_yvar-zsame_cohort1.rda"
dfs.path <- file.path("deconvo_method-paper", "outputs", folder.name, dfs.name)
dfs.new <- get(load(dfs.path))
# assign colnames
colnames(dfs.new) <- c("s.glial", "s.neuron", "s.train.variable.label", "s.train.variable.name")
# set numeric df for runs
s.col.index <- 1:2
for(c in s.col.index){dfs.new[,c] <- as.numeric(dfs.new[,c])}

#---------------------------
# get results for each s set
#---------------------------
# this is the chunk that makes the results df (CHECK CRUCIAL NOTES)
validation.sample.id.vector <- validation.sample.id.vector[validation.sample.id.vector %in% sce$Sample]
df.res.samples <- multigroup_bias_matched(validation.sample.id.vector[1:2], list.df.true, y.validate, dfs.new, sce)

# append s scale factor info
# assign dfs.new variable
dfs.new$s.value.string <- paste0(dfs.new$s.glial,";",dfs.new$s.neuron)
# get df.res.samples vector
s.value.vector <- paste0(df.res.samples$s.glial,";",df.res.samples$s.neuron)
# get new df.res variable
new.dfres.s.variable.name <- new.dfres.s.variable.label <- rep("NA", nrow(df.res.samples))
for(s.value.string in unique(s.value.vector)){
  message(s.value.string)
  # filter dfs
  filter.dfs <- dfs.new$s.value.string==s.value.string
  dfs.new.filtered <- dfs.new[filter.dfs,]
  # assign new values
  variable.value <- paste0(unique(dfs.new.filtered$variable.name), collapse = ";")
  variable.label <- paste0(unique(dfs.new.filtered$label), collapse = ";")
  # update new df.res variables
  filter.dfres.variable <- s.value.vector==s.value.string
  new.dfres.s.variable.name[filter.dfres.variable] <- variable.value
  new.dfres.s.variable.label[filter.dfres.variable] <- variable.label
}
df.res.samples$s.variable.name <- new.dfres.s.variable.name
df.res.samples$s.variable.label <- new.dfres.s.variable.label

for(value in s.value.vector){
  new.s.variable <- c()
}

df.res.samples$s.variable.label <- rep(dfs.new$label, each = ncol(y.validate))
df.res.samples$s.variable.name <- rep(dfs.new$variable.name, each = ncol(y.validate))

# append coldata from y.unadj (see MAE data)
y.unadj <- y.validate
cd.ydata <- colData(y.unadj)
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

# save/export
folder.name <- "13_soptimize_yvary-zvary_dlpfc-cohort1"
save.filename <- "df-sopt-result-validation_yvary-zvary_cohort1.rda"
save.path <- file.path("deconvo_method-paper", "outputs", folder.name, save.filename)
save(df.res.samples, file = save.path)
