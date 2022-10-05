#!/usr/bin/env R

# Plot mean ratios (MRs) and mean ratio summaries.

library(ggplot2)
library(gridExtra)

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
  theme(axis.text.x = element_blank(), legend.position = "none",
        axis.title.x = element_blank(), axis.title.y = element_blank()) +
  ggtitle("Violin plots")
ggjt2 <- ggplot(dfpi, aes(x = celltype.target, y = ratio, color = celltype.target)) + 
  theme_bw() + ylab("MR_k") + xlab("cell_type") + geom_jitter(size = 1, alpha = 0.5) +
  theme(axis.text.x = element_blank(), legend.position = "none",
        axis.title.x = element_blank(), axis.title.y = element_blank()) +
  ggtitle("Jitter points")
ggbp2 <- ggplot(dfpi, aes(x = celltype.target, y = ratio, fill = celltype.target)) + 
  theme_bw() + geom_boxplot() + ylab("MR_k") + xlab("cell_type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.position = "none") +
  ggtitle("Boxplots")

# save individual plots

# save composite plot
lm <- matrix(c(1,1,2,2,3,3,3),ncol=1) # adds extra space for xaxis labels

plot.comp.fname <- "mr-k_composite-vp-jt-bp.png"
plot.fpath <- file.path(save.dpath, plot.comp.fname)

png(plot.fpath, width = 4, height = 6, units = 'in', res = 400)

grid.arrange(ggvp2, ggjt2, ggbp2, ncol = 1, layout_matrix = lm, 
             left = paste0(paste0(rep(" ", 20),collapse=""), "MR_k"), 
             bottom = "cell_type")

dev.off()

# save composite plot

 







