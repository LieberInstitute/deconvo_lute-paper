#!/usr/bin/env R

# Author: Sean Maden
#
# Runs the main within-samples deconvolution trails, saving the results.

source("deconvo_method-paper/code/14_deconvolution-framework-trials/00_parameters.R")
sapply(libv, library, character.only = T)
lexperiment <- get(load(lexperiment.withinsample.path))
rse <- get(load(rse.k2markers.filepath))
image.table <- get(load(halo.output.path))
sce <- get(load(sce.markers.list.path))[["k2"]]
sce <- logNormCounts(sce, assay.type = "counts_adj")

# get matched datasets
# image
image.sample.id.vector <- standard_sample_id(image.table, "Position")
sce.sample.id.vector <- standard_sample_id(colData(sce), "Position")
rse.sample.id.vector <- standard_sample_id(colData(rse), "location")
# unique sample ids
unique.sample.id.sce <- unique(sce.sample.id.vector)
unique.sample.id.rse <- unique(rse.sample.id.vector)
unique.sample.id.image <- unique(image.sample.id.vector)
# overlapping sample ids
complete.sample.id.vector <- intersect(unique.sample.id.sce, 
                                       intersect(unique.sample.id.image, 
                                                 unique.sample.id.rse))

# perform deconvolution experiments
# get main z signature matrix
z.main <- signature_matrix_from_sce(sce)
s.main <- c("glial" = 3, "neuron" = 10)
# perform experiments
results.table <- do.call(rbind, lapply(complete.sample.id.vector, function(sample.id){
  lsample <- lexperiment[[sample.id]]
  z <- z.main; y <- lsample$y; s <- s.main
  do.call(rbind, lapply(method.vector, function(method){
    results.matrix <- do.call(rbind, lapply(seq(ncol(y)), function(index){
      yi <- y[,index,drop=FALSE]
      param.text <- paste0(method, "(y = yi, z = z, s = s)")
      new.parameters <- eval(parse(text = param.text))
      result <- deconvolution(new.parameters)
    }))
    results.matrix <- results.matrix %>% as.data.frame()
    # bind bulk experiment metadata
    results.matrix$bulk.id <- lsample$rse[["expt_condition"]]
    results.matrix$library.prep <- lsample$rse[["library_prep"]]
    results.matrix$library.type <- lsample$rse[["library_type"]]
    # bind proportions
    results.proportions <- apply(results.matrix[,c(1:2)], 1, function(row.index){
      row.index/sum(row.index)
    }) %>% t() %>% as.data.frame()
    colnames(results.proportions) <- paste0(colnames(results.proportions), 
                                            "_proportion")
    colnames(results.matrix)[1:2] <- paste0(colnames(results.matrix)[1:2], "_value")
    results.matrix <- cbind(results.matrix, results.proportions)
    results.matrix$sample.id <- sample.id
    results.matrix$method <- gsub("Param", "", method)
    results.matrix$glial.proportion.true <- lexperiment[[sample.id]]$p.proportion.k2["glial"]
    results.matrix$glial.count.true <- lexperiment[[sample.id]]$p.count.k2["glial"]
    results.matrix$neuron.proportion.true <- lexperiment[[sample.id]]$p.proportion.k2["neuron"]
    results.matrix$neuron.count.true <- lexperiment[[sample.id]]$p.count.k2["neuron"]
    return(results.matrix)
  }))
}))
# append errors, absolute errors
results.table$error.neuron <- results.table$neuron_proportion-
  results.table$neuron.proportion.true
results.table$absolute.error.neuron <- abs(results.table$error.neuron)
# append rmse
results.table$rmse <- apply(results.table[,c(6,7,10,12)], 1, function(row.index){
  sqrt(sum(row.index[1:2]-row.index[3:4]^2)/length(row.index[1:2]))
})
# save results table
save(results.table, file = within.samples.results.table.path)