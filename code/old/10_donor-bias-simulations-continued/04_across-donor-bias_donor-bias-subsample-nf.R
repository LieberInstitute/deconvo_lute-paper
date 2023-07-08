#!/usr/bin/env R

# Author: Sean Maden
#
# Run donor bias subsampling experiments.
#

source("deconvo_method-paper/code/10_donor-bias-simulations-continued/00_parameters.R")
sapply(libv, library, character.only = T)

# set up a new subsample experiment
sce <- get(load(sce.list.path))[[celltype.variable]]; rm(lscef)
proportions <- sce[[k.marker.variable]] %>% table() %>% prop.table()
lsub <- prepare_subsample_experiment(sce, scale.factor = scale.factor,
                                     iterations = iterations, groups.per.iteration = 3,
                                     method.vector = methodv, celltype.variable = celltype.variable,
                                     group.variable = group.variable, assay.name = assay.name,
                                     cell.proportions = proportions, count.minimum = count.min,
                                     seed.num = seed.num, which.save = c("sce", "tp", "ypb", "li"),
                                     save.fnstem = save.fnstem, base.directory.workflow.data = "data", verbose = TRUE)
setwd(rnf.workflow.acrossgroup.directory)
# save workflow table
workflow.table <- lsub$workflow.table
write.csv(workflow.table, file = workflow.table.path)
# save table iterations
index.vector <- seq(1, nrow(workflow.table), num.batch)
workflow.table.filename.ter <- paste0("workflow-table_",save.fnstem,".csv")
for(iteration.number in seq(length(index.vector))){
  index.start <- index.vector[iteration.number]
  index.end <- (index.start+num.batch-1)
  if(index.end > nrow(workflow.table)){index.end <- nrow(workflow.table)}
  workflow.tablei <- workflow.table[index.start:index.end,]
  workflow.tablei.fname.iter <- workflow.table.fnamei <- paste0("workflow-table-iter", 
                                        iter.num, "_", save.fnstem,".csv")
  workflow.tablei.fpath <- file.path(base.directory.workflow.data, workflow.tablei.fname.iter)
  write.csv(workflow.tablei, file = workflow.tablei.fpath, row.names = F)
  message("finished with iter ", iter.num)
}

# run the workflow
# navigate to:
# /deconvo_method-paper/outputs/10_donor-bias-simulations-continued/r-nf_deconvolution_donor-bias
#
# run:
# bash ./sh/r-nf_run_workflow.table-iter.sh