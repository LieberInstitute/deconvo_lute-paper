#!/usr/bin/env R

# Author: Sean Maden
#
# Test independent pseudobulk from multi-region brain data.

source("deconvo_method-paper/code/02_pseudobulk-outputs/00_parameters_script-set-14.R")
sapply(libv, library, character.only = T)
sce.mrb <- get(load(sce.mrb.path))

# get predictions



sce <- get(load(sce.markers.list.path))[["k2"]]

set.seed(0)

# prepare pseudobulk experiment
# filter sce.mrb
filter.markers <- rownames(sce.mrb) %in% rownames(sce)
sce.mrb <- sce.mrb[filter.markers,]
mrb.cell.type.vector <- sce.mrb[["cellType"]]
mrb.is.neuron <- grepl("Excit|Inhib", mrb.cell.type.vector)
mrb.is.glial <- grepl("Oligo|OPC|Micro|Astro", mrb.cell.type.vector)
colData(sce.mrb)[,"k2"] <- ifelse(mrb.is.neuron, "neuron", 
                                  ifelse(mrb.is.glial, "glial", "other"))
filter.cell.type <- sce.mrb[["k2"]] %in% c("neuron", "glial")
sce.mrb <- sce.mrb[,filter.cell.type]
# get logcounts
sce.mrb <- logNormCounts(sce.mrb, assay.type = "counts")
sce <- logNormCounts(sce, assay.type = "counts")
# cell size adjustments for pseudobulk
s.pb <- c("glial" = 2, "neuron" = 20)
# get pseudobulks for each donor
list.pb <- pseudobulk_from_sce(sce = sce.mrb, group.variable = "donor", s = s.pb, 
                               cell.type.variable = "k2", assay.name = "logcounts")
# get main signature matrix
z <- signature_matrix_from_sce(sce)

# experiment series 1 -- no cell size adjustment on z
# set cell size scales for z
s <- c("glial" = 1, "neuron" = 1)
# perform experiment
lresult.nnls.noadj <- run_pseudobulk_experiment(list.pb, s = s, z = z, method = "nnlsParam")
lresult.music.noadj <- run_pseudobulk_experiment(list.pb, s = s, z = z, method = "musicParam")
lresult.deconrnaseq.noadj <- run_pseudobulk_experiment(list.pb, s = s, z = z, method = "deconrnaseqParam")

# experiment series 1 -- with correct cell size adjustment on z
# set cell size scales for z
s <- s.pb
# perform experiment
lresult.nnls.adj <- run_pseudobulk_experiment(list.pb, s = s, z = z, method = "nnlsParam")
lresult.music.adj <- run_pseudobulk_experiment(list.pb, s = s, z = z, method = "musicParam")
lresult.deconrnaseq.adj <- run_pseudobulk_experiment(list.pb, s = s, z = z, method = "deconrnaseqParam")

# aggregate results table
donor.id.vector <- names(lresult.nnls.adj)
results.table <- do.call(rbind, lapply(donor.id.vector, function(donor.id){
  # results -- adj
  nnls.adj <- lresult.nnls.adj[[donor.id]]
  music.adj <- lresult.music.adj[[donor.id]]
  decon.adj <- lresult.deconrnaseq.adj[[donor.id]]
  # results -- noadj
  nnls.noadj <- lresult.nnls.noadj[[donor.id]]
  music.noadj <- lresult.music.noadj[[donor.id]]
  decon.noadj <- lresult.deconrnaseq.noadj[[donor.id]]
  
  # get row data
  row.nnls.adj <- c(nnls.adj$p.predicted, nnls.adj$bias, 
                    nnls.adj$rmse, "nnls", TRUE)
  row.music.adj <- c(music.adj$p.predicted, music.adj$bias, 
                     music.adj$rmse, "music", TRUE)
  row.decon.adj <- c(decon.adj$p.predicted, decon.adj$bias, 
                     decon.adj$rmse, "decon", TRUE)
  row.nnls.noadj <- c(nnls.noadj$p.predicted, nnls.noadj$bias, 
                      nnls.noadj$rmse, "nnls", FALSE)
  row.music.noadj <- c(music.noadj$p.predicted, music.noadj$bias, 
                       music.noadj$rmse, "music", FALSE)
  row.decon.noadj <- c(decon.noadj$p.predicted, decon.noadj$bias, 
                       decon.noadj$rmse, "decon", FALSE)
  
  # bind results table
  results <- rbind(row.nnls.adj, row.music.adj, row.decon.adj, 
                   row.nnls.noadj, row.music.noadj, row.decon.noadj) %>% as.data.frame()
  colnames(results) <- c("glial.prop.pred", "neuron.prop.pred", "glial.error",
                         "neuron.error", "rmse.k2", "method", "cell.size.adjustment")
  # append sample metadata
  results$sample.id <- donor.id
  results$neuron.true.prop <- list.pb[[donor.id]]$p[["neuron"]]
  results$glial.true.prop <- list.pb[[donor.id]]$p[["glial"]]
  results$neuron.true.count <- list.pb[[donor.id]]$p.counts[["neuron"]]
  results$glial.true.count <- list.pb[[donor.id]]$p.counts[["glial"]]
  return(results)
}))
results.table <- results.table %>% as.data.frame()
# format as data frame
for(index in c(1,2,3,4,5,9,10,11,12)){
  results.table[,index] <- results.table[,index] %>% as.numeric()}
# save results table
save(results.table, file = independent.pb.results.table.path)
