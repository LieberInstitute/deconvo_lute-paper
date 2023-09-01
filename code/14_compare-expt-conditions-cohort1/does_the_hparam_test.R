
#
#
#
#

#
#

libv <- c("ggplot2", "dplyr", "scuttle", "MultiAssayExperiment", "SingleCellExperiment")
sapply(libv, library, character.only = T)
knitr::opts_chunk$set(echo = TRUE)
folder.name <- "14_compare-expt-conditions-cohort1"

setwd("..")
setwd("..")
env.path <- file.path("outputs", folder.name, 
                      "02_only-rnascope-sizes_compare-experiment-conditions-cohort1.RData")
load(env.path)














# Gets search for single condition, single hyperparam set
expt.state <- "Nuc_polyA"
expt.var <- "expt_condition"
assay.name <- "logcounts"
group.name <- "sample.id"
yf.train <- y.train[,y.train[[expt.var]]==expt.state]
yf.validate <- y.validate[,y.validate[[expt.var]]==expt.state]
sample.id.vector <- c(unique(yf.train$sample.id)[1], unique(y.validate$sample.id)[1])
# search hyperparameters
sf.rn <- dfs.rn[dfs.rn$sample.id.dfs.rn == yf.train$sample.id[1],]
num.steps <- 35
s.step.validate <- 25
search.range.increment <- 8
s.search.min <- min(sf.rn[1:2])-search.range.increment
s.search.max <- max(sf.rn[1:2])+search.range.increment
s.search.step <- (s.search.max-s.search.min)/num.steps
dfsf.rn <- get_dfs(num.types = 2, min.size = s.search.min, 
                   max.size = s.search.max, size.step = s.search.step)
dfs.train <- dfsf.rn
colnames(dfs.train)[seq(2)] <- colnames(dfs.rn)[seq(2)]
dim(dfs.train)
# run experiment
list.expt <- get_soptimize_data_list(sce = sce, 
                                     sample.id.vector = sample.id.vector, 
                                     list.df.true = list.df.true, dfs = dfs.train,
                                     y.eset = y.unadj, y.train = y.train, 
                                     y.validate = y.validate, 
                                     assay.name = assay.name, 
                                     celltype.variable = k.variable.name, 
                                     group.name = group.name,
                                     matched.sce = TRUE)







###
# source(script.path) #  source latest script
###
crossval.unmatched.result <- 
  crossvalidate_soptimization(list.expt, s.step.validate = s.step.validate, 
                              validate.dfs = TRUE, plot.results = FALSE)
min(crossval.unmatched.result$df.res.train$error.glial.true.pred)
min(crossval.unmatched.result$df.res.validate$error.glial.true.pred)



