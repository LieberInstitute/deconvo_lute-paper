#!/usr/bin/env R

# Author: Sean Maden
#
# Run pseudobulk shuffle experiment varying Spseudobulk.
#

source(file.path("scripts","02_pseudobulk","00_param.R"))

# dependencies
knitr::opts_chunk$set(echo = TRUE)
libv <- c("snow", "dplyr", "parallel", "doParallel", "lute", "dplyr", "scuttle", "MultiAssayExperiment")
sapply(libv, library, character.only = T)

# 
#----------
# load data
#----------
# load mae (SEE CODE 01 OUTPUTS)
mae.filename <- "mae_allsamples.rda"
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
# df.rnascope.kdata
# df.rn.path <- file.path("outputs", "00_preprocess", "df-rnascope-info_cohort1.rda")
# df.rn <- get(load(df.rn.path))
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
sce <- sce[,sce[[sample.id.variable]] %in% sample.id.train]

#-----------------------------------
#
# experiment: high neuron proportion
#
#-----------------------------------
sample.id.iter <- sample.id.iter.high <- "Br3942_mid"
s.vector <- dfs.rn[dfs.rn$sample.id==sample.id.iter,]
s.vector <- c("glial" = as.numeric(s.vector[1]), "neuron" = as.numeric(s.vector[2]))
dfp.tall <- get_ypb_experiment_series(sce.k2, 
                                      sample.id.variable = sample.id.variable,
                                      celltype.variable = celltype.variable, 
                                      assay.name = assay.name,
                                      s.vector = s.vector,
                                      algorithm.name = algorithm.name, 
                                      return.dimensions = return.dimensions)
dfp.tall$s.sample.id <- sample.id.iter
dfp.tall$matched.id <- dfp.tall$s.sample.id==dfp.tall$sample.id
dfp.tall.high <- dfp.tall

#-----------------------------------
#
# experiment: low neuron proportion
#
#-----------------------------------
sample.id.iter <- sample.id.iter.low <- "Br2743_ant"
s.vector <- dfs.rn[dfs.rn$sample.id==sample.id.iter,]
s.vector <- c("glial" = as.numeric(s.vector[1]), "neuron" = as.numeric(s.vector[2]))
dfp.tall <- get_ypb_experiment_series(sce.k2, 
                                      sample.id.variable = sample.id.variable,
                                      celltype.variable = celltype.variable, 
                                      assay.name = assay.name,
                                      s.vector = s.vector,
                                      algorithm.name = algorithm.name, 
                                      return.dimensions = return.dimensions)
dfp.tall$s.sample.id <- sample.id.iter
dfp.tall$matched.id <- dfp.tall$s.sample.id==dfp.tall$sample.id
dfp.tall.low <- dfp.tall




# save

save.image(file = file.path("env", "03_shuffle", "00_fig3ab_script.RData"))

