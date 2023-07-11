#!/usr/bin/env R

#
# Testing relationship between cell sizes and total cells.
#

# aggregation basis code:
# do.call(rbind, lapply(sample.id.vector, function(sample.id){}))

source("deconvo_method-paper/code/05_cell-count-versus-cell-size/00_parameters.R")
sapply(libv, library, character.only = T)

# load mae 
mae.filepath <- here("deconvo_method-paper", "outputs", "01_prepare-datasets", "mae_final.rda")
mae <- get(load(mae.filepath))
# get mae metadata for subsetting
mae.cd <- colData(mae.final)
sample.id.vector <- mae.cd$sample.id[complete.cases(mae.final)]

# get table of cell sizes and cell counts
df.cellsize <- do.call(rbind, lapply(sample.id.vector, function(sample.id){
  mae.iter <- mae[,mae.cd$sample.id==sample.id,]
  
  # gets cell size
  img.data <- mae.iter[[3]] %>% assays() %>% unlist() %>% t() %>% as.data.frame()
  img.data$cell_type <- mae.iter[[3]]$cell_type
  df.img <- aggregate(img.data[,"Nucleus_Area"], list("cell_type" = img.data$cell_type), FUN = "median")
  colnames(df.img)[2] <- "median_nucleus_area"
  
  # gets cell amount by type, from sce
  stat <- sumCountsAcrossCells(mae.iter[[2]], ids = colData(mae.iter[[2]])[,"cellType_broad_hc"])
  stat.med <- colMedians(assays(stat)[["sum"]])
  df.sn <- data.frame(median_counts_snrnaseq = stat.med, cell_type = colnames(stat))
  # gets total cell amount by type, from sce
  stat2 <- sumCountsAcrossCells(mae.iter[[2]], ids = rep("all", ncol(mae.iter[[2]])))
  stat.med2 <- colMedians(assays(stat2)[["sum"]])
  df.sn.all <- data.frame(median_counts_snrnaseq = stat.med2, cell_type = colnames(stat2))
  df.sn.final <- rbind(df.sn, df.sn.all)
 
  # harmonized image and sn data
  df.harmonized <- harmonize_celltype_tables_1to1(df.sn.final, df.img)
  
  df.harmonized <- append_k_columns(df.harmonized)
  df.harmonized$sample.id <- sample.id
  df.harmonized
}))

#-----------------------------------------
# plots -- snrnaseq vs. image nucleus area
#-----------------------------------------
dfp <- df.cellsize[,c(1,2,4,6:9)]
dfp$median_counts_snrnaseq <- as.numeric(dfp$median_counts_snrnaseq)
dfp$median_nucleus_area <- as.numeric(dfp$median_nucleus_area)
dfp$cell_type <- as.character(dfp$cell_type)
dfp$sample.id <- as.character(dfp$sample.id)
dfp <- dfp[!is.na(dfp$median_nucleus_area),]

# facet by sample id
ggplot(dfp, aes(x = median_counts_snrnaseq, y = median_nucleus_area, color = cell_type)) +
  geom_point() + facet_wrap(~sample.id) + geom_abline(slope = 1, intercept = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(dfp, aes(x = median_counts_snrnaseq, y = median_nucleus_area, shape = cell_type)) +
  geom_point(size=4,alpha=0.5) + facet_wrap(~sample.id) + geom_abline(slope = 1, intercept = 0) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(dfp, aes(x = median_counts_snrnaseq, y = median_nucleus_area, shape = k2)) +
  geom_point(size=4,alpha=0.5) + facet_wrap(~sample.id) + geom_abline(slope = 1, intercept = 0) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(dfp, aes(x = median_counts_snrnaseq, y = median_nucleus_area, shape = k3)) +
  geom_point(size=4,alpha=0.5) + facet_wrap(~sample.id) + geom_abline(slope = 1, intercept = 0) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(dfp, aes(x = median_counts_snrnaseq, y = median_nucleus_area, shape = k4)) +
  geom_point(size=4,alpha=0.5) + facet_wrap(~sample.id) + geom_abline(slope = 1, intercept = 0) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(dfp, aes(x = median_counts_snrnaseq, y = median_nucleus_area, color = sample.id)) +
  geom_point() + facet_wrap(~cell_type) + geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# glial only
ggplot(dfp[dfp$k2=="glial",], aes(x = median_counts_snrnaseq, y = median_nucleus_area, color = sample.id)) +
  geom_point() + facet_wrap(~cell_type) + geom_abline(slope = 1, intercept = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(dfp[dfp$k2=="glial",], aes(x = median_counts_snrnaseq, y = median_nucleus_area)) +
  geom_point(alpha = 0.5) + facet_wrap(~cell_type) + geom_abline(slope = 1, intercept = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme_bw()

# k2 labels
ggplot(dfp, aes(x = median_counts_snrnaseq, y = median_nucleus_area)) + theme_bw() +
  geom_point(alpha = 0.5) + facet_wrap(~k2) + geom_abline(slope = 1, intercept = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# k3 labels
ggplot(dfp, aes(x = median_counts_snrnaseq, y = median_nucleus_area)) + theme_bw() +
  geom_point(alpha = 0.5) + facet_wrap(~k3) + geom_abline(slope = 1, intercept = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# k4 labels
ggplot(dfp, aes(x = median_counts_snrnaseq, y = median_nucleus_area)) + theme_bw() +
  geom_point(alpha = 0.5) + facet_wrap(~k4) + geom_abline(slope = 1, intercept = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

