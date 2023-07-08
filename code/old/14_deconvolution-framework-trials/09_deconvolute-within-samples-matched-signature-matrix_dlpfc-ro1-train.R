
# perform deconvolution
# uses matched pseudobulk and signature matrices within samples
# includes null cell sizes adjustment

source("deconvo_method-paper/code/14_deconvolution-framework-trials/00_parameters_script-set-14.R")
sapply(libv, library, character.only = T)
lexperiment <- get(load(lexperiment.withinsample.path))

# image.table <- get(load(halo.output.path))
sce <- get(load(sce.markers.list.path))[["k2"]]
sce <- logNormCounts(sce, assay.type = "counts")
# complete.sample.id.vector <- get(load(complete.sample.id.vector.path))

# perform deconvolution experiments
# get main z signature matrix
z.main <- z <- signature_matrix_from_sce(sce)
# perform experiments
markers.per.type <- 20
cell.type.variable <- "k2"
assay.name <- "counts"
sample.id.vector <- names(lexperiment)
results.list <- lapply(sample.id.vector, function(sample.id){
  lsample <- lexperiment[[sample.id]]
  y.sample <- lsample$y
  sizes.list <- lsample$list.sizes
  #z.sample <- signature_matrix_from_sce(lsample$sce, 
  #                               cell.type.variable = cell.type.variable,
  #                               summary.method = "mean", 
  #                               assay.name = assay.name)
  z.sample <- lute(lsample$sce, 
                   y = y.sample,
                   markers.per.type = markers.per.type, 
                   assay.name = assay.name,
                   celltype.variable = cell.type.variable,
                   typemarker.algorithm = "meanratio",
                   deconvolution.algorithm = NULL)
  
  sample.results <- results_series_table(y = y.sample, z = z.sample, 
                                         method.vector = method.vector,
                                         sizes.list = sizes.list)
  sample.results$bulk.label <- colnames(y.sample)
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