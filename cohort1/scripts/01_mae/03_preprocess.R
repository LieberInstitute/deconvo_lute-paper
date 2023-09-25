#!/usr/bin/env R

#
# Preprocess MultiAssayExperiment
#

min.neuron.proportion <- 0.2
max.nucleus.area <- 78

#-----
# load
#-----

mae.in.path <- "./outputs/01_mae/mae_allsamples_append.rda"
mae <- get(load(mae.in.path))
cd <- colData(mae)

#----------------------------------------------
# 1. filter sample IDs on cell type proportions
#----------------------------------------------
list.true.proportions <- metadata(
  mae[["snrnaseq.k2.all"]])[["list.df.true.k2"]]

filter.condition <- sapply(
  list.true.proportions, function(item){
    item[["neuron"]] > min.neuron.proportion})

identical(list.true.proportions[filter.condition], list.true.proportions)

sample.id.exclude <- names(list.true.proportions[!filter.condition])

filter.mae <- !colData(mae)$sample.id %in% sample.id.exclude
mae <- mae[,filter.mae,]
dim(mae)

#------------------------------------------
# 2. filter RNAscope/HALO cells on max area
#------------------------------------------
rn <- mae[["cell.sizes"]]

#-----
# save
#-----

mae.out.path <- "./outputs/01_mae/mae_analysis.rda"
mae <- get(load(mae.path))

