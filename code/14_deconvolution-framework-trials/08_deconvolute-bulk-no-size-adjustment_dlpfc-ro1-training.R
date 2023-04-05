
# Get the unadjusted deconvolution results for figure 1b, using real bulk samples.
# Matched signature matrix and bulk data from DLPFC RO1 training.

source("deconvo_method-paper/code/14_deconvolution-framework-trials/00_parameters_script-set-14.R")
sapply(libv, library, character.only = T)
lexperiment <- get(load(lexperiment.withinsample.path))
image.table <- get(load(halo.output.path))
sce <- get(load(sce.markers.list.path))[["k2"]]
sce <- logNormCounts(sce, assay.type = "counts")

# load bulk data
bulk.name <- "rse_k2-marker-expression_ro1-dlpfc.rda"
bulk.path <- here("deconvo_method-paper", "outputs", "11_bulk-expression-analysis")
bulk.path <- here(bulk.path, bulk.name)
bulk.rse <- get(load(bulk.path))
bulk.counts <- assays(bulk.rse)[["counts"]]
# set sample id colnames
bulk.sample.id.vector <- colnames(bulk.counts)
bulk.sample.id.vector <- paste0("Br", gsub(".*Br", "", bulk.sample.id.vector))
colnames(bulk.counts) <- bulk.sample.id.vector

# set cell sizes list to be equal (not applied)
lsize.no.adjustment <- list("none" = c("glial" = 1, "neuron" = 1))

# get deconvolution results series
results.list <- lapply(sample.id.vector, function(sample.id){
  lsample <- lexperiment[[sample.id]]
  
  
  y <- bulk.expression[,sample.id] # lsample$y
  sizes.list <- lsample$list.sizes
  sample.results <- results_series_table(y = y, z = z, 
                                         method.vector = method.vector,
                                         sizes.list = lsize.no.adjustment)
  
  # sample.results$bulk.label <- colnames(y)
  sample.results$sample.id <- sample.id
  sample.results$y.total.expression <- lsample$y.total.expression
  sample.results$glial.proportion.true <- lsample$p.proportion.k2["glial"]
  sample.results$glial.count.true <- lsample$p.count.k2["glial"]
  sample.results$neuron.proportion.true <- lsample$p.proportion.k2["neuron"]
  sample.results$neuron.count.true <- lsample$p.count.k2["neuron"]
  return(sample.results)
})
results.table <- do.call(rbind, results.list) %>% as.data.frame()