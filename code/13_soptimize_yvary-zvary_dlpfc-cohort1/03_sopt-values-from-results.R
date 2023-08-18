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
load.filename <- "df-sopt-result_yvary-zvary_cohort1.rda"
load.path <- file.path("deconvo_method-paper", "outputs", folder.name, load.filename)
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


#------------------
# save results list
#------------------
# get results list to save
lr <- list(results = df.res, minima.by.label = df.min, sopt.values.by.label = list.sopt)

# saves minima with results
lr.filename <- "list-sopt-values_results-yvar-zvary_cohort1.rda"
results.list.path <- file.path("deconvo_method-paper", "outputs", folder.name, lr.filename)
save(lr, file = results.list.path)

#----------------------------------------------------
# save median sopt.train by condition -- UNSUPERVISED
#----------------------------------------------------
# load
df.min <- lr$minima.by.label
# source function: dfs_byvariable()
folder.name <- "13_soptimize_yvary-zvary_dlpfc-cohort1"
script.path <- file.path("deconvo_method-paper", "code", folder.name, "00_parameters.R")
source(script.path)
# assign categories
df.min$compartment_library <- df.min$cell.compartment
df.min$library.preparation <- gsub(".*_", "", df.min$cell.compartment)
df.min$cell.compartment <- gsub("_.*", "", df.min$cell.compartment)
# get new dfs as medians by experiment group category
variable.vector <- c("anatomic.region", "cell.compartment", "compartment_library", "library.preparation", "sample.id")
dfs.new <- dfs_byvariable(df.min, variable.vector)

# save 
dfs.name <- "dfs-medians-bygroup-training_yvar-zsame_cohort1.rda"
dfs.path <- file.path("deconvo_method-paper", "outputs", folder.name, dfs.name)
save(dfs.new, file = dfs.path)

#-------------------------------------------------------
# save median sopt.train by condition -- SEMI-SUPERVISED
#-------------------------------------------------------
# load
df.min <- lr$minima.by.label
# source function: dfs_byvariable()
folder.name <- "13_soptimize_yvary-zvary_dlpfc-cohort1"
script.path <- file.path("deconvo_method-paper", "code", folder.name, "00_parameters.R")
source(script.path)
# assign categories
df.min$compartment_library <- df.min$cell.compartment
df.min$library.preparation <- gsub(".*_", "", df.min$cell.compartment)
df.min$cell.compartment <- gsub("_.*", "", df.min$cell.compartment)
# get new dfs as medians by experiment group category
variable.vector <- c("anatomic.region", "cell.compartment", "compartment_library", "library.preparation", "sample.id")
dfs.new <- dfs_byvariable(df.min, variable.vector)

# save 
dfs.name <- "dfs-medians-bygroup-training_SUPERVISED-FILTER_yvar-zsame_cohort1.rda"
dfs.path <- file.path("deconvo_method-paper", "outputs", folder.name, dfs.name)
save(dfs.new, file = dfs.path)

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



