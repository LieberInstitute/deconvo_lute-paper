#!/usr/bin/env R

# Author: Sean Maden
#
# Test independent pseudobulk from multi-region brain data.

source("deconvo_method-paper/code/14_deconvolution-framework-trials/00_parameters.R")
sapply(libv, library, character.only = T)
sce.mrb <- get(load(sce.mrb.path))
sce <- get(load(sce.markers.list.path))[["k2"]]

# prepare deconvolution experiment
s <- c("glial" = 5, "neuron" = 10)
# filter sce.mrb
filter.markers <- rownames(sce.mrb) %in% rownames(sce)
sce.mrb <- sce.mrb[filter.markers,]
mrb.cell.type.vector <- sce.mrb[["cellType"]]
mrb.is.neuron <- grepl("Excit|Inhib", mrb.cell.type.vector)
mrb.is.glial <- grepl("Oligo|OPC|Micro|Astro", mrb.cell.type.vector)
colData(sce.mrb)[,"k2"] <- ifelse(mrb.is.neuron, "neuron", 
                                  ifelse(mrb.is.glial, "glial", "other"))
filter.cell.type <- sce.mrb[["k2"]] %in% c("neuron", "glial")
sce.mrb <- sce.mrb[,filter.cell.type]
# get logcounts
sce.mrb <- logNormCounts(sce.mrb)
sce <- logNormCounts(sce)
# get pseudobulks for each donor
list.pb <- pseudobulk_from_sce(sce = sce.mrb, group.variable = "donor", 
                               s = s, cell.type.variable = "k2", 
                               assay.name = "logcounts")
# get main signature matrix
z <- signature_matrix_from_sce(sce)

# perform experiment
lresult.nnls <- run_experiment(list.pb, method = "nnlsParam")
lresult.music <- run_experiment(list.pb, method = "musicParam")