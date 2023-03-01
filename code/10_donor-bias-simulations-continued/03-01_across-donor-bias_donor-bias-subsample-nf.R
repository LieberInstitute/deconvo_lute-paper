#!/usr/bin/env R

# Author: Sean Maden
#
# Run donor bias subsampling experiments.
#
# Note: run after script before 03-02 to use the same across-sample bulk 
# reference across the two experiments.
#

libv <- c("lute", "SummarizedExperiment", "SingleCellExperiment", 
          "ggplot2", "gridExtra")
sapply(libv, library, character.only = TRUE)

#------------------
# experiment params
#------------------
# get load path
code.dname <- "09_manuscript"
proj.dname <- "deconvo_method-paper"
load.dpath <- file.path(proj.dname, "outputs", code.dname)
# get save path
code.dname <- "10_donor-bias-simulations-continued"
proj.dname <- "deconvo_method-paper"
save.dpath <- file.path(proj.dname, "outputs", code.dname)
# get params for experiment
seed.num <- 0
celltype.variable <- "k2"
proj.handle <- "ro1-dlpfc"
save.fnstem <- paste0("inter-sample_", proj.handle)
group.variable <- "Sample"
assay.name <- "counts_adj"
methodv <- c("nnls", "music", "epic", "deconrnaseq")
iterations <- 1000
fraction.cells <- 25
num.sample.iter <- 3
scale.factor <- c("glial" = 3, "neuron" = 10)
rnf.dname <- "r-nf_deconvolution_donor-bias"
base.path <- "data"
base.path <- file.path(save.dpath, rnf.dname, base.path)

#----------------------------------
# set up a new subsample experiment
#----------------------------------
# load data
fname <- paste0("list-scef_markers-k2-k3-k4_",proj.handle,".rda")
fpath <- file.path(load.dpath, fname)
lscef <- get(load(fpath))
sce <- lscef[[celltype.variable]]
rm(lscef)

# save new experiment data
lsub <- prepare_subsample_experiment(sce,
                                       scale.factor = scale.factor,
                                       iterations = iterations,
                                       groups.per.iteration = 3,
                                       method.vector = methodv,
                                       celltype.variable = celltype.variable,
                                       group.variable = group.variable,
                                       assay.name = assay.name,
                                       fraction.cells = fraction.cells,
                                       seed.num = seed.num,
                                       which.save = c("sce", "tp", "ypb", "li"),
                                       save.fnstem = save.fnstem,
                                       base.path = "data",
                                       verbose = TRUE)

#---------------------
# manage workflow runs
#---------------------
# change wd
# setwd(file.path(save.dpath, rnf.dname))

# get main starting workflow table
wt <- lsub$wt
# save full table
proj.handle <- "ro1-dlpfc"
save.fnstem <- paste0("inter-sample_", proj.handle)
wt.fnamei <- paste0("workflow-table-all_",save.fnstem,".csv")
wt.fpath <- file.path(save.dpath, rnf.dname, "data", wt.fnamei)
write.csv(wt, file = wt.fpath)

# save table iterations
num.batch <- 200
indexv <- seq(1, nrow(wt), num.batch)
wt.fnamei <- paste0("workflow-table_",save.fnstem,".csv")
for(iter.num in seq(length(indexv))){
  indexv.iter <- indexv[iter.num]
  indexv.iter <- indexv.iter:(indexv.iter+num.batch-1)
  wti <- wt[indexv.iter,]
  wti.fname.iter <- wt.fnamei <- paste0("workflow-table-iter", 
                                       iter.num, "_", save.fnstem,".csv")
  wti.fpath <- file.path(save.dpath, rnf.dname, "data", wti.fname.iter)
  write.csv(wti, file = wti.fpath, row.names = F)
  message("finished with iter ", iter.num)
}

#-----------------------
# analyze results table
#-----------------------
results.filt <- "results_table_.*"
data.dpath <- file.path(save.dpath, rnf.dname, "data")
lfv <- list.files(data.dpath)
rt.fname <- lfv[grepl(results.filt, lfv)]
rt.fpath <- file.path(save.dpath, rnf.dname, "data", rt.fname)
rt <- read.csv(rt.fpath)



