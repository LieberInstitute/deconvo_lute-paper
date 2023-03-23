#!/usr/bin/env R

# Author: Sean Maden
#
# Run subsampling experiments within batch groups.
# 

source("deconvo_method-paper/code/10_donor-bias-simulations-continued/00_parameters.R")
sapply(libv, library, character.only = T)

# set up a new subsample experiment
sce <- get(load(sce.list.path))[[celltype.variable]]; rm(lscef)
# filter on largest donor
group.table <- sce[[group.variable]] %>% table() %>%  as.data.frame()
filter.group.table <- dft[,2]==max(dft[,2])
largest.donor.id <- dft[filter.group.table, 1]
sce.filter <- sce[[group.variable]]==largest.donor.id
sce.filtered <- sce[,sce.filter]
# prepare workflow table
proportions <- scef[[celltype.variable.name]] %>% table() %>% prop.table()
lsub <- prepare_subsample_experiment(sce.filtered,
                                     scale.factor = scale.factor,
                                     iterations = iterations,
                                     groups.per.iteration = 1,
                                     method.vector = methodv,
                                     celltype.variable = celltype.variable,
                                     group.variable = group.variable,
                                     assay.name = assay.name,
                                     cell.proportions = proportions,
                                     count.minimum.acrossgroupimum = count.minimum.acrossgroup,
                                     seed.num = seed.num,
                                     which.save = c("sce", "tp", "ypb", "li"),
                                     save.fnstem = save.fnstem,
                                     base.path = "data",
                                     verbose = TRUE)

# manage workflow runs
setwd(workflow.path.withingroup)
workflow.table <- lsub$workflow.table
write.csv(workflow.table, file = workflow.table.path.withingroup)

# save table iterations
index.vector <- seq(1, nrow(workflow.table), num.batch.withingroup)
workflow.table.fnamei <- paste0("workflow-table_",save.fnstem,".csv")
for(iteration.number in seq(length(index.vector))){
  index.start <- index.vector[iteration.number]
  index.end <- (index.start+num.batch.withingroup-1)
  if(index.end > nrow(workflow.table)){index.end <- nrow(workflow.table)}
  workflow.tablei <- workflow.table[index.start:index.end,]
  workflow.tablei.fname.iter <- workflow.table.fnamei <- paste0("workflow-table-iter", iteration.number, "_", save.fnstem,".csv")
  workflow.tablei.fpath <- file.path(base.path, workflow.tablei.fname.iter)
  write.csv(workflow.tablei, file = workflow.tablei.fpath, row.names = F)
  message("finished with iter ", iteration.number)
}

# run the workflow
# navigate to:
# /deconvo_method-paper/outputs/10_donor-bias-simulations-continued/r-nf_deconvolution_donor-bias
#
# run:
# bash ./sh/r-nf_run_workflow.table-iter.sh