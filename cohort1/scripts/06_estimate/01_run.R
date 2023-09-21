#!/usr/bin/env R

# Author: Sean Maden
#
# Get S estimates from bulk validation samples, as a result of multi-step S optimization.
#

libv <- c("snow", "dplyr", "parallel", "doParallel", "lute", "dplyr", "MultiAssayExperiment")
sapply(libv, library, character.only = T)

# sets variables
script.path <- file.path("scripts", "06_estimate", "00_param.R")
source(script.path)
source("./scripts/06_estimate/00_dataset_summaries.R")
source("./scripts/06_estimate/00_deconvo_plots.R")

#----------
# load data
#----------
# load mae (SEE CODE 01 OUTPUTS)
new.mae.filename <- "mae_allsamples.rda"
mae.final.filepath <- file.path("outputs", "01_mae", new.mae.filename)
mae <- get(load(mae.final.filepath))

# true proportions from sce data
sce.all <- mae[["snrnaseq.k2.all"]]
sample.id.vector <- unique(sce.all[["Sample"]])
list.df.true <- lapply(sample.id.vector, function(sample.id){
  k.table <- table(sce.all[,sce.all[["Sample"]]==sample.id][["k2"]])
  k.prop <- as.data.frame(t(as.matrix(prop.table(k.table))))
  rownames(k.prop) <- "true_proportion"
  return(k.prop)
})
names(list.df.true) <- sample.id.vector

#-------------------------------
# subset validation and training
#-------------------------------
train.sample.id <- get(load("./outputs/00_preprocess/list_snrnaseq_sampleid.rda"))[["train"]]
cd.mae <- colData(mae)
filter.mae.train <- cd.mae$sample.id %in% train.sample.id
mae.train <- mae[,filter.mae.train]
mae.validate <- mae[,!filter.mae.train]
rm(mae)

#-----------
# experiment
#-----------
# run experiment
seq.steps.train <- seq(5, 100, 5)
list.cv <- lapply(seq.steps.train, function(train.steps){
  get_crossvalidation_results(
    mae.train, mae.validate, num.steps.train = train.steps,
    num.steps.validate = 10)
})

# get plot data
dfp <- do.call(rbind, lapply(list.cv, function(cv.result){
  c(min(cv.result$validate.result$df.res$error.neuron), cv.result$num.steps.train)
}))
dfp <- as.data.frame(dfp)
colnames(dfp) <- c("min.error.neuron.validate", "num.steps.train")

#-----
# save
#-----
# files

list.cv.validate <- lapply(list.cv, function(item){
  list(results.validate = item$validate.result, num.steps.train = item$num.steps.train)
})
save(list.cv.validate, file = "./outputs/06_estimate/results_cv_validate.rda")

save(list.cv, file = "./outputs/06_estimate/results_cv.rda")
save(dfp, file = "./outputs/06_estimate/results_dfp.rda")

# env
rm(mae.final)
rm(mae.validate)
rm(mae.train)
save.image("./outputs/06_estimate/01_run_script.RData")

#-------------
# plot results
#-------------
ggplot(dfp, aes(x = num.steps.train, y = min.error.neuron.validate)) + geom_point()

# plot
dfp <- data.frame(min.neuron.error.val = c(min(cv10$validate.result$df.res$error.neuron),
                                           min(cv20$validate.result$df.res$error.neuron)),
                  steps.train = c(10, 20))



cv5 <- get_crossvalidation_results(
  mae.train, mae.validate, num.steps.train = 5)

cv10 <- get_crossvalidation_results(
  mae.train, mae.validate, num.steps.train = 10)

cv15 <- get_crossvalidation_results(
  mae.train, mae.validate, num.steps.train = 15)

cv20 <- res20 <- get_crossvalidation_results(
  mae.train, mae.validate, num.steps.train = 20)


# plot
dfp <- data.frame(min.neuron.error.val = c(min(cv10$validate.result$df.res$error.neuron),
                                           min(cv20$validate.result$df.res$error.neuron)),
                  steps.train = c(10, 20))
ggplot(dfp, aes(x = steps.train, y = min.neuron.error.val)) + geom_point()


cv50 <- get_crossvalidation_results(
  mae.train, mae.validate, num.steps.train = 50)

cv100 <- get_crossvalidation_results(
  mae.train, mae.validate, num.steps.train = 100)






save(list.dfres.train, file = "./")




# set params (SEE PROJECT NOTES)
assay.name <- "counts"
celltype.variable <- "k2"
sample.id.variable <- "Sample"
sample.id.vector <- colData(mae.train)$sample.id
filter.sample.id <- sample.id.vector %in% mae.train[["bulk.rnaseq"]][[y.group.name]]
sample.id.vector <- sample.id.vector[filter.sample.id]
y.group.name <- 'batch.id2'
dfs.iter <- dfs.series()
nrow(dfs.iter)
# iterate on training samples
list.res <- lapply(sample.id.vector, function(sample.id){
  filter.mae.train <- colData(mae.train)$sample.id==sample.id
  mae.train.iter <- mae.train[,filter.mae.train,]
  message("Num. tests: ", nrow(dfs.iter))
  # seq
  y.iter <- mae.train.iter[["bulk.rnaseq"]]
  sce.iter <- mae.train.iter[["snrnaseq.k2.all"]]
  # results
  multigroup_bias_matched(sample.id, list.df.true, y.iter, 
                          y.group.name = y.group.name,
                          dfs.iter, sce.iter, 
                          assay.name = assay.name)
})
df.res.train <- do.call(rbind, lapply(list.res, function(item){item}))
df.res.train$crossvalidation <- "train"
summary(df.res.train$bias.neuron.true.pred)
# prepare and plot results
df.res.train <- dfres_postprocess(df.res.train)
list.plots.dfres.train <- get_dfres_plots(df.res)

# save results
save(df.res.train, file = "./outputs/06_estimate/train_result.rda")
save(list.plots.dfres.train, file = "./outputs/06_estimate/train_plot.rda")

# get dfs
dfs.iter <- dfs.series(seq(1, 20, 0.5))
nrow(dfs.iter)




#--------------------
# s optimize validate
#--------------------
# get dfs from train min.error coordinates
min.error.neuron <- min(df.res.train$error.neuron)
min.error.neuron # 7.89636e-05
filter.res <- df.res.train$error.neuron == min.error.neuron
s.train.neuron <- df.res.train[filter.res,]$s.neuron
s.train.glial <- df.res.train[filter.res,]$s.glial
s.vector.validate <- c(s.train.neuron, s.train.glial)
s.validate.min <- min(s.vector.validate)
s.validate.max <- max(s.vector.validate)
s.num.steps <- 400
s.validate.increment <- (s.validate.max-s.validate.min)/s.num.steps
s.validate.seq <- seq(s.validate.min, s.validate.max, s.validate.increment)
dfs.validate <- dfs.series(s.validate.seq)

# iterate on training samples
sample.id.vector <- unique(
  intersect(
    mae.validate[["bulk.rnaseq"]][["batch.id2"]], 
    mae.validate[["snrnaseq.k2.all"]][["Sample"]]))
list.res.validate <- lapply(sample.id.vector, function(sample.id){
  filter.mae <- colData(mae.validate)$sample.id==sample.id
  mae.iter <- mae.validate[,filter.mae,]
  message("Num. tests: ", nrow(dfs.iter))
  # seq
  y.iter <- mae.iter[["bulk.rnaseq"]]
  sce.iter <- mae.iter[["snrnaseq.k2.all"]]
  # results
  multigroup_bias_matched(sample.id, list.df.true, y.iter, 
                          y.group.name = y.group.name,
                          dfs.iter, sce.iter, 
                          assay.name = assay.name)
})
df.res.validate <- do.call(rbind, lapply(list.res.validate, function(item){item}))
df.res.validate$crossvalidation <- "validate"
summary(df.res.validate$error.neuron)

for(sample.id in unique(df.res.validate$sample.id)){
  filt <- df.res.validate$sample.id==sample.id
  df.res.validate[filt,]$neuron.true <- list.df.true[[sample.id]]$neuron
  df.res.validate[filt,]$glial.true <- list.df.true[[sample.id]]$glial
}
df.res.validate$bias.neuron.true.pred <- df.res.validate$neuron.true-df.res.validate$neuron.pred.nnls
df.res.validate$bias.glial.true.pred <- df.res.validate$glial.true-df.res.validate$glial.pred.nnls

# prepare and plot results
df.res <- dfres_postprocess(df.res.validate)
list.plots.dfres <- get_dfres_plots(df.res)

# plot
ggplot(df.res, aes(x = sample.id, y = error.neuron)) + 
  geom_violin(draw_quantiles=0.5)

#----------------------
# get s summary results
#----------------------
# rnascope

# s optimize

# marker library size













# plot error by sample id
ggplot(df.res.train, aes(x = sample.id, y = error.neuron)) + 
  geom_violin(draw_quantiles = 0.5)

ggplot(df.res.train[filter.train,], aes(x = sample.id, y = error.neuron)) + 
  geom_violin(draw_quantiles = 0.5)

barplot(table(df.res.train[filter.train,]$sample.id))

ggplot(table(df.res.train[filter.train,]$sample.id), 
       aes(x = sample.id)) + 
  geom_barplot()

ggplot(df.res.train[filter.train,], 
       aes(x = s.glial, y = s.neuron)) + 
  geom_point()

# plot
ggplot(df.res.train[filter.train,], aes(x = 1, y = s.neuron)) + 
  geom_violin(draw_quantiles = 0.5)
ggplot(df.res.train[filter.train,], aes(x = 1, y = s.glial)) + 
  geom_violin(draw_quantiles = 0.5)

df.res.train$s.neuron


train.min.errors <- lapply(sample.id.vector, function(sample.id){
  dff <- df.res.train[df.res.train$sample.id==sample.id,]
  min.error.neuron <- min(abs(dff$bias.neuron.true.pred))
  message("min. neuron error: ", min.err.neuron)
  if(min.error.neuron <= max.error){
    filter.dff <- abs(dff$bias.neuron.true.pred)==min.error.neuron
    return(c(median(dff.min$s.neuron), 
             median(dff.min$s.glial)))
  }
})



# save
df.res.train <- df.res.samples
list.plots.dfres.train <- list.plots.dfp1
save(df.res.train, file = "./outputs/06_estimate/train_result.rda")
save(list.plots.dfres.train, file = "./outputs/06_estimate/train_plot.rda")

#-----
# save
#-----
