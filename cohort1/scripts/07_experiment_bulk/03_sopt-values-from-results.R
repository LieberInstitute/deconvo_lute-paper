#!/usr/bin/env R

# Author: Sean Maden
#
# Gets the optimal s cell size factor values from run results.
#

library(ggplot2)

#-----
# load
#-----
load.filename <- "df-sopt-result_yvary-zsame_cohort1.rda"
load.path <- file.path("deconvo_method-paper", "outputs", 
                       "12_soptimize_yvary-zsame_dlpfc-cohort1", load.filename)
df.res <- get(load(load.path))

# get unique labels
labels.vector <- unique(df.res$sample.label)

#---------------------
# get optimal s values
#---------------------
# gets the minima
df.min <- do.call(rbind, lapply(labels.vector, function(label.iter){
  df.iter <- df.res[df.res$sample.label==label.iter,]
  df.min.iter <- df.iter[df.iter$error.neuron == min(df.iter$error.neuron),,drop=F]
  df.min.iter[1,]
}))
df.min$min.error.neuron <- df.min$error.neuron

# get list of optimal s values by labels
list.sopt <- lapply(labels.vector, function(label.iter){
  df.iter <- df.min[df.min$sample.label==label.iter,]
	c("glial" = df.iter[,"s.glial"], "neuron" = df.iter[,"s.neuron"])
})
names(list.sopt) <- labels.vector

#---------------
# plot summaries
#---------------

dfp <- rbind(data.frame(s.train = df.min$s.glial), data.frame(s.train = df.min$s.neuron))
dfp$cell.type <- c(rep("glial", nrow(df.min)), rep("neuron", nrow(df.min)))

# violin
ggplot(dfp, aes(x = cell.type, y = s.train)) + geom_violin(draw_quantiles = 0.5)

# boxplot
ggplot(dfp, aes(x = cell.type, y = s.train)) + 
  geom_jitter(alpha = 0.5) + geom_boxplot(color = "cyan", alpha = 0)

# scatter
ggplot(df.min, aes(x = s.glial, y = s.neuron)) + geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0)
ggplot(df.min, aes(x = s.glial, y = s.neuron)) + geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0) + facet_wrap(~cell.compartment)

#------------------
# save results list
#------------------
# get results list to save
lr <- list(results = df.res, minima.by.label = df.min, sopt.values.by.label = list.sopt)

# saves minima with results
lr.filename <- "list-sopt-values_results-yvar_cohort1.rda"
results.list.path <- file.path("deconvo_method-paper", "outputs", 
                               "12_soptimize_yvary-zsame_dlpfc-cohort1", lr.filename)
save(lr, file = results.list.path)







