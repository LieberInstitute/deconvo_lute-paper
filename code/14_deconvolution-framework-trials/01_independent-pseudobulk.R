#!/usr/bin/env R

# Author: Sean Maden
#
# Test independent pseudobulk from multi-region brain data.

source("deconvo_method-paper/code/14_deconvolution-framework-trials/00_parameters.R")
sapply(libv, library, character.only = T)
sce.mrb <- get(load(sce.mrb.path))
sce <- get(load(sce.markers.list.path))[["k2"]]

# prepare deconvolution experiment
s <- c("glial" = 5, "neuron" = 10)
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
sce.mrb <- logNormCounts(sce.mrb)
sce <- logNormCounts(sce)
# get pseudobulks for each donor
list.pb <- pseudobulk_from_sce(sce = sce.mrb, group.variable = "donor", 
                               s = s, cell.type.variable = "k2", 
                               assay.name = "logcounts")
# get main signature matrix
z <- signature_matrix_from_sce(sce)

# perform experiment
lresult.nnls <- run_pseudobulk_experiment(list.pb, method = "nnlsParam")
lresult.music <- run_pseudobulk_experiment(list.pb, method = "musicParam")
lresult.deconrnaseq <- run_pseudobulk_experiment(list.pb, method = "deconrnaseqParam")

# aggregate results table
donor.id.vector <- names(lresult.nnls)
results.table <- do.call(rbind, lapply(donor.id.vector, function(donor.id){
  nnls <- lresult.nnls[[donor.id]]
  music <- lresult.music[[donor.id]]
  decon <- lresult.deconrnaseq[[donor.id]]
  # get row data
  row.nnls <- c(nnls$p.predicted, nnls$bias, nnls$rmse, "nnls")
  row.music <- c(music$p.predicted, music$bias, music$rmse, "music")
  row.decon <- c(decon$p.predicted, decon$bias, decon$rmse, "deconrnaseq")
  # bind results table
  results <- rbind(row.nnls, row.music, row.decon) %>% as.data.frame()
  colnames(results) <- c("glial.prop.pred",
                         "neuron.prop.pred",
                         "glial.error",
                         "neuron.error",
                         "rmse.k2",
                         "method")
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
for(index in c(1,2,3,4,5,8,9,10,11)){
  results.table[,index] <- as.numeric(results.table[,index])}
# save results table
save(results.table, file = independent.pb.results.table.path)

# plot true vs. pred proportions, neuron
new.plot <- ggplot(results.table, aes(x = neuron.true.prop, 
                                      y = neuron.prop.pred, 
                                      color = method)) +
  geom_point(alpha = 0.6) + geom_abline(slope = 1, intercept = 0, color = "black") +
  geom_vline(xintercept = 0.5, alpha = 0.5) + xlab("True") + ylab("Predicted") +
  ggtitle("Cell type: neuron") + xlim(0, 1) + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

new.plot + facet_wrap(~sample.id)
jpeg(pb.scatterplot.proportions.bysample.colmethod.path, width = 5, 
     height = 4, units = "in", res = 400)
new.plot + facet_wrap(~sample.id)

# plot absolute errors, neuron
results.table$neuron.abs.error <- abs(results.table$neuron.error)
new.plot <- ggplot(results.table, 
                   aes(x = method, y = neuron.abs.error, fill = method)) +
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  xlab("Deconvolution algorithm") +
  ylab("Absolute error") +
  ggtitle("Cell type: neuron")
jpeg(pb.barplot.abserror.bysample.colmethod.path, width = 5, height = 4, 
     units = "in", res = 400)
new.plot + facet_wrap(~sample.id)
dev.off()

# plot rmse, k2
new.plot <- ggplot(results.table, aes(x = method, y = rmse.k2, fill = method)) +
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  xlab("Deconvolution algorithm") +
  ylab("RMSE") +
  ggtitle("Cell type: neuron")
jpeg(pb.barplot.rmse.bysample.colmethod.path, width = 5, height = 4, 
     units = "in", res = 400)
new.plot + facet_wrap(~sample.id)
dev.off()
