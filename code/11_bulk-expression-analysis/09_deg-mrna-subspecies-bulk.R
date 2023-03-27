#!/usr/bin/env R

# Author: Sean Maden
#
# Compare expression of mRNA sub-species (histone, mitochondrial, etc.)

# load data
source("deconvo_method-paper/code/11_bulk-expression-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
rse.background <- get(load(rse.gene.filter.filepath))
cd <- colData(rse.background)
rd <- rowData(rse.background)
table(rd$gene_type) # check types
# lncRNA        Mt_rRNA        Mt_tRNA protein_coding           rRNA 
# 17748              2             22          19988             47

# get rse subsets
# histone genes
# get histone gene expression
histone.filter <- grepl(histone.genes.pattern, rd$Symbol)
rse.histone <- rse.background[histone.filter,]
dim(rse.histone)
# [1]  97 113
# get mtrna
mito.filter <- grepl(mito.genes.pattern, rd$Symbol)
rse.mito <- rse.background[mito.filter,]
dim(rse.mito)
# [1]  37 113
# ribo rna filter
ribo.filter <- rd$gene_type == "rRNA"
rse.rrna <- rse.background[ribo.filter,]
dim(rse.rrna)
# [1]  47 113
# lincRNA filter
lncrna.filter <- rd$gene_type == "lncRNA"
rse.linc <- rse.background[lncrna.filter,]
dim(rse.linc)
# [1] 17748   113
# expression set list, with gene counts
lrse <- list(Histone = rse.histone, Mitochondrial = rse.mito, 
             Linc = rse.linc, Ribosomal = rse.rrna)
gene.counts <- unlist(lapply(lrse, function(rse){nrow(rse)}))
names(lrse) <- paste0(names(lrse), " (", gene.counts," genes)")

# get plot data
# means
data.mean <- do.call(rbind, lapply(seq(length(lrse)), function(index){
  rse <- lrse[[index]]; cd <- colData(rse); expression <- assays(rse)[[assay.name]]
  plot.data <- group_jitter("expt_condition", cd, expression, "mean")$dfp
  plot.data$type <- names(lrse)[index]
  return(plot.data)
}))
mean.plot <- ggplot(data.mean, aes(x = group, y = value)) + theme_bw() +
  geom_jitter(alpha = 0.8) + geom_boxplot(alpha = 0, color = "cyan") +
  facet_wrap(~type) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab(paste0("Mean expression (log10 ", assay.name, ")")) +
  xlab("Experiment group") + scale_y_log10()
# variances
data.variance <- do.call(rbind, lapply(seq(length(lrse)), function(index){
  rse <- lrse[[index]]; cd <- colData(rse); expression <- assays(rse)[[assay.name]]
  plot.data <- group_jitter("expt_condition", cd, expression, "variance")$dfp
  plot.data$type <- names(lrse)[index]
  return(plot.data)
}))
variance.plot <- ggplot(data.variance, aes(x = group, y = value)) + theme_bw() +
  geom_jitter(alpha = 0.8) + geom_boxplot(alpha = 0, color = "cyan") +
  facet_wrap(~type) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab(paste0("Expression variance (log10 ", assay.name, ")")) +
  xlab("Experiment group") + scale_y_log10()
# total counts
data.total <- do.call(rbind, lapply(seq(length(lrse)), function(index){
  rse <- lrse[[index]]; cd <- colData(rse); expression <- assays(rse)[[assay.name]]
  plot.data <- group_jitter("expt_condition", cd, expression, "total.counts")$dfp
  plot.data$type <- names(lrse)[index]
  return(plot.data)
}))
total.plot <- ggplot(data.total, aes(x = group, y = value)) + theme_bw() +
  geom_jitter(alpha = 0.8) + geom_boxplot(alpha = 0, color = "cyan") +
  facet_wrap(~type) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab(paste0("Total expression (log10 ", assay.name, ")")) +
  xlab("Experiment group") + scale_y_log10()

# save new plots
plot.width <- 6
plot.height <- 5
plot.res <- 400
# save
jpeg(bulk.mean.rna.types.jpg.path, width = plot.width, 
     height = plot.height, units = "in", res = plot.res)
mean.plot
dev.off()
# save
jpeg(bulk.variance.rna.types.jpg.path, width = plot.width, 
     height = plot.height, units = "in", res = plot.res)
variance.plot
dev.off()
# save
jpeg(bulk.total.rna.types.jpg.path, width = plot.width, 
     height = plot.height, units = "in", res = plot.res)
total.plot
dev.off()

# composite plots
# mean expression
# jpeg(mean.expression.rna.types.jpg.name, width = 10, height = 5, units = "in", res = 400)
