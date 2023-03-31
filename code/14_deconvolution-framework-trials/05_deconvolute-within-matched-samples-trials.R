#!/usr/bin/env R

# Author: Sean Maden
#
# Runs the main within-samples deconvolution trails, saving the results.

source("deconvo_method-paper/code/14_deconvolution-framework-trials/00_parameters.R")
sapply(libv, library, character.only = T)
lexperiment <- get(load(lexperiment.withinsample.path))
image.table <- get(load(halo.output.path))
sce <- get(load(sce.markers.list.path))[["k2"]]
sce <- logNormCounts(sce, assay.type = "counts")
# complete.sample.id.vector <- get(load(complete.sample.id.vector.path))

# perform deconvolution experiments
# get main z signature matrix
z.main <- z <- signature_matrix_from_sce(sce)
# perform experiments
sample.id.vector <- names(lexperiment)
results.list <- lapply(sample.id.vector, function(sample.id){
  lsample <- lexperiment[[sample.id]]
  y <- lsample$y
  sizes.list <- lsample$list.sizes
  sample.results <- results_series_table(y = y, z = z, 
                                         method.vector = method.vector,
                                         sizes.list = sizes.list)
  sample.results$bulk.label <- colnames(y)
  sample.results$sample.id <- sample.id
  sample.results$y.total.expression <- lsample$y.total.expression
  sample.results$glial.proportion.true <- lsample$p.proportion.k2["glial"]
  sample.results$glial.count.true <- lsample$p.count.k2["glial"]
  sample.results$neuron.proportion.true <- lsample$p.proportion.k2["neuron"]
  sample.results$neuron.count.true <- lsample$p.count.k2["neuron"]
  return(sample.results)
})
results.table <- do.call(rbind, results.list) %>% as.data.frame()
# append errors, absolute errors
results.table$error.neuron <- results.table$neuron.proportion.predicted-
  results.table$neuron.proportion.true
results.table$absolute.error.neuron <- abs(results.table$error.neuron)
# append rmse
# save results table
save(results.table, file = within.samples.results.table.path)