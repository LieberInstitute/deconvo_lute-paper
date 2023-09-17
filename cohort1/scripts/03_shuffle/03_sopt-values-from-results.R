#!/usr/bin/env R

# Author: Sean Maden
#
# Gets the optimal s cell size factor values from TRAINING run results.
#

libv <- c("ggplot2")
sapply(libv, library, character.only = T)

#-----
# load
#-----
folder.name <- "13_soptimize_yvary-zvary_dlpfc-cohort1"
# load train results
load.filename <- "df-sopt-result_yvary-zvary_cohort1.rda"
load.path <- file.path("deconvo_method-paper", "outputs", folder.name, load.filename)
df.res <- get(load(load.path))

# source function: dfs_byvariable()
script.path <- file.path("deconvo_method-paper", "code", folder.name, "00_parameters.R")
source(script.path)

#---------------------
# get optimal s values
#---------------------
# gets data.frame of error minima from sopt.training results
labels.vector <- unique(df.res$sample.label)
df.min <- do.call(rbind, lapply(labels.vector, function(label.iter){
  df.iter <- df.res[df.res$sample.label==label.iter,]
  df.min.iter <- df.iter[df.iter$error.neuron == min(df.iter$error.neuron),,drop=F]
  df.min.iter[1,]
}))
df.min$min.error.neuron <- df.min$error.neuron

# differentiate df.min data for sopt summaries
df.min.unsupervised <- df.min
filter.supervision <- df.min$s.glial < df.min$s.neuron
df.min.supervised <- df.min[filter.supervision,]

# get sopt median summaries
# note: get new dfs as medians by experiment group category
variable.vector <- c("anatomic.region", "cell.compartment", 
                     "cell.compartment.library.type", 
                     "library.type", "sample.id")
dfs.sopt.unsupervised <- dfs_byvariable(df.min.unsupervised, variable.vector)
dfs.sopt.supervised <- dfs_byvariable(df.min.supervised, variable.vector)
dfs.sopt.unsupervised$train.type <- "unsupervised"
dfs.sopt.supervised$train.type <- "supervised"
dfs.sopt.train <- rbind(dfs.sopt.unsupervised, dfs.sopt.supervised)

# save
dfs.name <- "dfs-medians-bygroup-training_yvary-zvary_cohort1.rda"
dfs.path <- file.path("deconvo_method-paper", "outputs", folder.name, dfs.name)
save(dfs.sopt.train, file = dfs.path)

#---------------
# plot summaries
#---------------
df.min <- df.min.supervised
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
