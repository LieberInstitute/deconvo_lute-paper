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
df.min <- list.sopt.train$minima.by.label
# assign categories
df.min$compartment_library <- df.min$cell.compartment
df.min$library.preparation <- gsub(".*_", "", df.min$cell.compartment)
df.min$cell.compartment <- gsub("_.*", "", df.min$cell.compartment)

# get new dfs as medians by experiment group category
variable.vector <- c("anatomic.region", "cell.compartment", "compartment_library", 
                     "library.preparation", "sample.id")
dfs.new <- dfs_byvariable(df.min, variable.vector)

# save 
dfs.name <- "dfs-medians-bygroup-training_yvar-zsame_cohort1.rda"
dfs.path <- file.path("deconvo_method-paper", "outputs", "12_soptimize_yvary-zsame_dlpfc-cohort1", dfs.name)
save(dfs.new, file = dfs.path)

# set numeric df for runs
dfs <- dfs.new[,c("s.glial", "s.neuron")]
for(c in seq(ncol(dfs))){dfs[,c] <- as.numeric(dfs[,c])}

#---------------------------
# get results for each s set
#---------------------------
# this is the chunk that makes the results df (CHECK CRUCIAL NOTES)
validation.sample.id.vector <- validation.sample.id.vector[validation.sample.id.vector %in% sce$Sample]
df.res.samples <- multigroup_bias_matched(validation.sample.id.vector, list.df.true, y.validate, dfs, sce)

# append coldata from y.unadj (see MAE data)
y.unadj <- y.validate
df.res.samples$sample.labels <- rep(colnames(y.unadj), nrow(dfs))
df.res.samples$sample.id <- rep(gsub("_.*", "", y.unadj$batch.id), nrow(dfs))
df.res.samples$cell.compartment <- rep(y.unadj$library_prep, nrow(dfs))
df.res.samples$anatomic.region <- rep(y.unadj$location, nrow(dfs))
df.res.samples$library.type <- rep(y.unadj$library_type, nrow(dfs))
df.res.samples$compartment_library <- rep(y.unadj$expt_condition, nrow(dfs))
df.res.samples$sample.id.brnum <- rep(y.unadj$batch.id2, nrow(dfs))
# append dfs info
df.res.samples$dfs.condition.label <- rep(dfs.new$label, ncol(y.unadj))
df.res.samples$dfs.condition.variable.name <- rep(dfs.new$variable.name, ncol(y.unadj))

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
save.filename <- "df-sopt-result-validation_yvary-zsame_cohort1.rda"
save.path <- file.path("deconvo_method-paper", "outputs", 
                       "12_soptimize_yvary-zsame_dlpfc-cohort1", save.filename)
save(df.res.samples, file = save.path)

#--------------
# example plots
#--------------
# compartment_library
condition_comparison_boxplots("compartment_library", "Bulk_RiboZeroGold", df.res.samples)
condition_comparison_boxplots("compartment_library", "Bulk_polyA", df.res.samples)
condition_comparison_boxplots("compartment_library", "Cyto_RiboZeroGold", df.res.samples)
condition_comparison_boxplots("compartment_library", "Cyto_polyA", df.res.samples)
condition_comparison_boxplots("compartment_library", "Nuc_RiboZeroGold", df.res.samples)
condition_comparison_boxplots("compartment_library", "Nuc_polyA", df.res.samples)

# library.preparation
df.res.samples$library.preparation <- df.res.samples$library.type
condition_comparison_boxplots("library.preparation", "polyA", df.res.samples)
condition_comparison_boxplots("library.preparation", "RiboZeroGold", df.res.samples)

# cell.compartment
condition_comparison_boxplots("cell.compartment", "Nuc", df.res.samples)
condition_comparison_boxplots("cell.compartment", "Cyto", df.res.samples)
condition_comparison_boxplots("cell.compartment", "Bulk", df.res.samples)

# anatomic.region
condition_comparison_boxplots("anatomic.region", "Ant", df.res.samples)
condition_comparison_boxplots("anatomic.region", "Mid", df.res.samples)
condition_comparison_boxplots("anatomic.region", "Post", df.res.samples)

# sample id
condition_comparison_boxplots("sample.id", "Br2720", df.res.samples)
condition_comparison_boxplots("sample.id", "Br6423", df.res.samples)
condition_comparison_boxplots("sample.id", "Br6432", df.res.samples)
