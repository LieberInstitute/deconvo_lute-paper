#!/usr/bin/env R

#
# Testing relationship between cell sizes and total cells.
#

# aggregation basis code:
# do.call(rbind, lapply(sample.id.vector, function(sample.id){}))

source("deconvo_method-paper/code/05_cell-count-versus-cell-size/00_parameters.R")
sapply(libv, library, character.only = T)

# load snrnaseq data
list.sce.markers <- get(load(sce.markers.list.path))
sce <- list.sce.markers$k2

#-------------------------------
# get cell sizes and counts data
#-------------------------------

sample.id.vector <- unique(sce$Sample)
assay.name <- "counts"

# get table of cell sizes and cell counts
df.cellsize <- do.call(rbind, lapply(sample.id.vector, function(sample.id){
  sce.iter <- sce[,sce$Sample == sample.id]
  # gets cell size -- this is the total expression by marker set
  expr.total.counts <- colSums(assays(sce.iter)[["counts"]])
  df.count <- data.frame(total.counts = expr.total.counts,
                         cell.type = sce.iter$cellType_broad_hc)
  df.count.agg <- aggregate(df.count$total.counts, by = list(df.count$cell.type), FUN = "median")
  colnames(df.count.agg) <- c("cell_type", "median_total_expression")
  # gets cell amount by type, from sce
  stat <- sumCountsAcrossCells(sce.iter, ids = sce.iter$cellType_broad_hc)
  # stat <- sumCountsAcrossCells(mae.iter[[2]], ids = colData(mae.iter[[2]])[,"cellType_broad_hc"])
  stat.med <- colMedians(assays(stat)[["sum"]])
  #df.sn <- data.frame(median_counts_snrnaseq = stat.med, cell_type = colnames(stat))
  df.sn.final <- data.frame(median_counts_snrnaseq = stat.med, cell_type = colnames(stat))
  # gets **total** cell amount by type, from sce
  #stat2 <- sumCountsAcrossCells(sce.iter, ids = rep("all", ncol(sce.iter)))
  #stat.med2 <- colMedians(assays(stat2)[["sum"]])
  #df.sn.all <- data.frame(median_counts_snrnaseq = stat.med2, cell_type = colnames(stat2))
  #df.sn.final <- rbind(df.sn, df.sn.all)
  
  # harmonized image and sn data
  #df.harmonized <- harmonize_celltype_tables_1to1(df.sn.final, df.img)
  #df.harmonized <- append_k_columns(df.harmonized)
  #df.harmonized$sample.id <- sample.id
  df.harmonized <- cbind(df.count.agg, df.sn.final)[,c(1:3)]
  df.harmonized$sample.id <- sample.id
  df.harmonized
}))
df.cellsize <- append_k_columns(df.cellsize)
df.cellsize <- df.cellsize[,c(1,2,3,4,6,7,8)]

#-----------------------------------------
# plots -- snrnaseq vs. image nucleus area
#-----------------------------------------
dfp <- df.cellsize
#dfp$median_counts_snrnaseq <- as.numeric(dfp$median_counts_snrnaseq)
#dfp$median_nucleus_area <- as.numeric(dfp$median_nucleus_area)
dfp$cell_type <- as.character(dfp$cell_type)
dfp$sample.id <- as.character(dfp$sample.id)
#dfp <- dfp[!is.na(dfp$median_nucleus_area),]

# facet by sample id
ggplot(dfp, aes(x = median_counts_snrnaseq, y = median_total_expression, color = cell_type)) +
  geom_point() + facet_wrap(~sample.id) + geom_abline(slope = 1, intercept = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(dfp, aes(x = median_counts_snrnaseq, y = median_total_expression, shape = cell_type)) +
  geom_point(size=4,alpha=0.5) + facet_wrap(~sample.id) + geom_abline(slope = 1, intercept = 0) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(dfp, aes(x = median_counts_snrnaseq, y = median_total_expression, shape = k2)) +
  geom_point(size=4,alpha=0.5) + facet_wrap(~sample.id) + 
  geom_abline(slope = 1, intercept = 0) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(dfp, aes(x = median_counts_snrnaseq, y = median_total_expression, shape = k3)) +
  geom_point(size=4,alpha=0.5) + facet_wrap(~sample.id) + 
  geom_abline(slope = 1, intercept = 0) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(dfp, aes(x = median_counts_snrnaseq, y = median_total_expression, shape = k4)) +
  geom_point(size=4,alpha=0.5) + facet_wrap(~sample.id) + 
  geom_abline(slope = 1, intercept = 0) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(dfp, aes(x = median_counts_snrnaseq, y = median_total_expression, color = sample.id)) +
  theme_bw() + geom_point() + facet_wrap(~cell_type) + geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
