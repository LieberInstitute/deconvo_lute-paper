#!/usr/bin/env R

#
# Get cell type proportions
#

libv <- c("dplyr", "ggplot2")
sapply(libv, library, character.only = TRUE)

source("./scripts/07_summary/00_size.R")

# load
mae <- get(load("./outputs/01_mae/mae_allsamples.rda"))
sn <- mae[["snrnaseq.k2.all"]]
sample.id.variable <- unique(sn[["Sample"]])
mae <- mae[,colData(mae)$sample.id %in% sample.id.variable,]
# sopt cell size factors
sopt <- get(load("./outputs/06_estimate/train_result.rda"))


#----------------
# get cell sizes
#----------------
list.dfp <- list_dfp_wide_tall_size(mae)

# get factors from optimizations
df.sopt <- do.call(rbind, lapply(unique(sopt$sample.id), function(sample.id){
  sf <- sopt[sopt$sample.id==sample.id,]
  c(sample.id, mean(sf$s.glial), mean(sf$s.neuron))
}))

# save env
rm(mae)
save.image("./env/07_summary/01_sizes_script.RData")

#---------------------------------------------
# compare size estimates within sample sources
#---------------------------------------------

