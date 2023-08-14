#!/usr/bin/env R

# Author: Sean Maden
#
# Gets the optimal s cell size factor values from run results.
#

#-----
# load
#-----
load.filename <- "df-sopt-result_yvary-zsame_cohort1.rda"
load.path <- file.path("outputs", "12_soptimize_yvary-zsame_dlpfc-cohort1", load.filename)
df.res <- get(load(load.path))
# get unique labels
labels.vector <- unique(df.res$sample.labels)


#---------------------
# get optimal s values
#---------------------
# gets the minima
df.min <- do.call(rbind, lapply(labels.vector, function(label.iter){
  df.iter <- df.res[df.res$sample.labels==label.iter,]
  df.iter[df.iter$error.neuron == min(df.iter$error.neuron),,drop=F]
}))
df.min$min.error.neuron <- df.min$error.neuron

# get list of optimal s values by labels
list.sopt <- lapply(labels.vector, function(label.iter){
	c("glial" = df.iter["glial"], "neuron" = df.iter["neuron"])
})
names(list.sopt) <- labels.vector

#------------------
# save results list
#------------------
# get results list to save
lr <- list(results = df.res, minima.by.label = df.min, sopt.values.by.label = list.sopt)
# saves minima with results
save(lr, lr.path)
