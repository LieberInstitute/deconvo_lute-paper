#!/usr/bin/env R

# Author: Sean Maden
#
# Run pseudobulk shuffle experiment varying 
#
#
#
#
#

new.env.name <- "03_figfg_script.RData"

source(file.path("scripts","02_pseudobulk","00_param.R"))
source(file.path("scripts","03_shuffle","00_param.R"))
source(file.path("scripts","03_shuffle","00_param_pseudobulk.R"))

# dependencies
knitr::opts_chunk$set(echo = TRUE)
libv <- c("snow", "dplyr", "parallel", "doParallel", "lute", "dplyr", "scuttle", "MultiAssayExperiment")
sapply(libv, library, character.only = T)

#----------
# load data
#----------
# load mae (SEE CODE 01 OUTPUTS)
mae.filename <- "mae_analysis_append.rda"
mae.path <- file.path("outputs", "01_mae", mae.filename)
mae <- get(load(mae.path))

folder.name <- "03_shuffle" 
celltype.variable <- "k2"
sample.id.variable <- "Sample"
assay.name <- "logcounts"
sce <- sce.k2 <- mae[["snrnaseq.k2.all"]]
algorithm.name <- "nnls"
return.dimensions <- "tall"

# load rnascope cell sizes
sample.id.variable.rn <- "sample.id"
df.rn <- mae[["cell.sizes"]]
dfs.rn <- t(df.rn) %>% as.data.frame()
dfs.rn <- dfs.rn[dfs.rn$k.label==celltype.variable,]
dfs.rn <- data.frame(s.glial = dfs.rn[dfs.rn$cell_type=="glial",]$cell_size,
                     s.neuron = dfs.rn[dfs.rn$cell_type=="neuron",]$cell_size,
                     sample.id.column = dfs.rn[dfs.rn$cell_type=="neuron",]$sample_id)
colnames(dfs.rn)[3] <- sample.id.variable.rn
for(c in seq(2)){dfs.rn[,c] <- as.numeric(dfs.rn[,c])}

# filter training samples
sample.id.train.path <- file.path("outputs", "00_preprocess", "list_snrnaseq_sampleid.rda")
sample.id.train <- get(load(sample.id.train.path))[["train"]]
sample.id.train <- sample.id.train[sample.id.train %in% dfs.rn$sample.id]
length(sample.id.train)

# subset on rnascope samples
sce <- sce[,sce[[sample.id.variable]] %in% sample.id.train]
sce.k2 <- sce.k2[,sce.k2$Sample %in% sample.id.train]

#-----------------------------------
#
# experiment: high neuron proportion
#
#-----------------------------------

zpb.sample.id <- zpb.sample.id.high <- "Br3942_mid"
sce.foriter <- sce.k2
sce.deconvo <- sce.k2[,sce.k2$Sample==zpb.sample.id]
dfp.tall <- get_ypb_experiment_series_shuffle_zpb(sce.foriter, 
                                                  sce.deconvo,
                                                  sample.id.variable = sample.id.variable,
                                                  celltype.variable = celltype.variable,
                                                  assay.name = assay.name,
                                                  algorithm.name = algorithm.name,
                                                  return.dimensions = return.dimensions)

dfp.tall$z.deconvo.sample.id <- zpb.sample.id
dfp.tall$matched.z <- dfp.tall$sample.id == dfp.tall$z.deconvo.sample.id
dfp.tall.high <- dfp.tall

#-----------------------------------
#
# experiment: low neuron proportion
#
#-----------------------------------

zpb.sample.id <- zpb.sample.id.low <- "Br2743_ant"
sce.foriter <- sce.k2
sce.deconvo <- sce.k2[,sce.k2$Sample==zpb.sample.id]
dfp.tall <- get_ypb_experiment_series_shuffle_zpb(sce.foriter, 
                                                  sce.deconvo,
                                                  sample.id.variable = sample.id.variable,
                                                  celltype.variable = celltype.variable,
                                                  assay.name = assay.name,
                                                  algorithm.name = algorithm.name,
                                                  return.dimensions = return.dimensions)

dfp.tall$z.deconvo.sample.id <- zpb.sample.id
dfp.tall$matched.z <- dfp.tall$sample.id == dfp.tall$z.deconvo.sample.id
dfp.tall.low <- dfp.tall


# save

save.image(file = file.path("env", "03_shuffle", new.env.name))
