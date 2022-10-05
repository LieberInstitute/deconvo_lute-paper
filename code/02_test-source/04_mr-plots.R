#!/usr/bin/env R

# Plot mean ratios (MRs) and mean ratio summaries.

library(ggplot2)

#------
# paths
#------
proj.dpath <- "deconvo_method-paper"
save.dpath <- file.path(proj.dpath, "outputs/02_test-source")
lz.fname <- "lz_mr_dlpfc-ro1.rda"

#-----
# load
#-----
lz <- get(load(file.path(save.dpath, lz.fname)))

#-------------
# violin plots
#-------------
dfpi <- lz$top.marker.data
dfpi <- as.data.frame(dfpi[,c(1,2,4,6,7)])

# mean ratio by cell type
# plots for individual figures
ggvp1 <- ggplot(dfpi, aes(x = celltype.target, y = ratio, fill = celltype.target)) + 
  theme_bw() + geom_violin(draw_quantiles = 0.5) + ylab("MR_k") + xlab("cell_type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggjt1 <- ggplot(dfpi, aes(x = celltype.target, y = ratio, color = celltype.target)) + 
  theme_bw() + ylab("MR_k") + xlab("cell_type") + geom_jitter(size = 1, alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggbp1 <- ggplot(dfpi, aes(x = celltype.target, y = ratio, fill = celltype.target)) + 
  theme_bw() + geom_boxplot() + ylab("MR_k") + xlab("cell_type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
# plots for composite figures
ggvp2 <- ggplot(dfpi, aes(x = celltype.target, y = ratio, fill = celltype.target)) + 
  theme_bw() + geom_violin(draw_quantiles = 0.5) + ylab("MR_k") + xlab("cell_type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggjt2 <- ggplot(dfpi, aes(x = celltype.target, y = ratio, color = celltype.target)) + 
  theme_bw() + ylab("MR_k") + xlab("cell_type") + geom_jitter(size = 1, alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggbp2 <- ggplot(dfpi, aes(x = celltype.target, y = ratio, fill = celltype.target)) + 
  theme_bw() + geom_boxplot() + ylab("MR_k") + xlab("cell_type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

# save individual plots


# save composite plot

 







