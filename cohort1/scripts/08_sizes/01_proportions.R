#!/usr/bin/env R

#
# Get cell type sizes
#

libv <- c("dplyr", "ggplot2")
sapply(libv, library, character.only = TRUE)

source("./scripts/07_summary/00_prop.R")

# load
mae <- get(load("./outputs/01_mae/mae_allsamples.rda"))
sn <- mae[["snrnaseq.k2.all"]]
sample.id.variable <- unique(sn[["Sample"]])
mae <- mae[,colData(mae)$sample.id %in% sample.id.variable,]

# params
min.proportion.thresh <- 0.2

#----------------
# get proportions
#----------------
list.dfp <- list_dfp_wide_tall(mae)

#----------------------------
# save min proportion filters
#----------------------------
df.wide <- list.dfp$dfp.wide
sn.id.filter.vector <- unique(df.wide[df.wide$sn.neuron < min.proportion.thresh,]$sample.id)
rn.id.filter.vector <- unique(df.wide[df.wide$rn.neuron < min.proportion.thresh,]$sample.id)
md.str <- paste0("sample.id vectors to exclude based on absolute proportions threshold of ", 
                 min.proportion.thresh)
list.proportion.filter <- list(rnascope = rn.id.filter.vector,
                               snrnaseq = sn.id.filter.vector,
                               min.proportion.thresh = min.proportion.thresh,
                               metadata = md.str)

#-----
# save
#-----
# save
save(list.proportion.filter, file = "./outputs/07_summary/list_proportions_filter.rda")

# save env
save.image("./env/07_summary/01_proportions_script.RData")
