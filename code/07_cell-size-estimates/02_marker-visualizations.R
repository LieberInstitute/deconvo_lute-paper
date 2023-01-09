#!/usr/bin/env R

# Author: Sean Maden
#
# Visualize k=2 markers
#
#

libv <- c("lute", "SingleCellExperiment", "SummarizedExperiment",
          "ComplexHeatmap", "ggplot2", "ggrepel")
sapply(libv, library, character.only = T)

#----------
# load data 
#----------
# se marker data, k2
sef.fname <- "sef_mr-markers_k2_20-per-k_dlpfc-ro1.rda"
sef.dpath <- file.path("deconvo_method-paper", "outputs", 
                       "05_marker-gene-annotations")
sef <- get(load(file.path(sef.dpath, sef.fname)))

# plot output paths
plot.dpath <- file.path("deconvo_method-paper", "outputs", "07_cell-size-estimates")

#---------------------
# get types expression
#---------------------
# get signature matrix, z
sef[["donor"]] <- sef[["BrNum"]]
setf <- set_from_sce(sef, groupvar = "donor", method = "mean",
                     assayname = "logcounts")
lct <- assays(setf)$logcounts

#------------------------------------
# visualizations -- marker expression
#------------------------------------
plot.fname.stem <- "markers-k2_n20-perk"
# heatmap of marker logcounts
plot.fname <- paste0("chm_",plot.fname.stem,".pdf")
pdf(file.path(plot.dpath, plot.fname), width = 4, height = 8)
Heatmap(assays(setf)$logcounts, name = "Logcounts", show_column_dend = F)
dev.off()

# make dfp
lct <- assays(setf)$logcounts
lct$marker <- rownames(lct)
dfp <- rbind(data.frame(mean = lct$Neuron, marker = lct$marker),
             data.frame(mean = lct$`Non-neuron`, marker = lct$marker))
dfp$celltype <- c(rep("Neuron", nrow(lct)), rep("Non-neuron", nrow(lct)))
dfp$mean <- round(dfp$mean, 2)

# tile plot
ggplot(dfp, aes(x = celltype, y = marker, label = mean, fill = mean)) + 
  geom_tile() + geom_text() + theme_bw()

# violin plot
ggplot(dfp, aes(x = celltype, y = mean, fill = celltype)) + 
  geom_violin(draw_quantiles = 0.5) + theme_bw()

# violin with scatter
plot.fname <- paste0("ggvp_",plot.fname.stem,".pdf")
pdf(file.path(plot.dpath, plot.fname), width = 2.5, height = 3)
ggplot(dfp, aes(x = celltype, y = mean, fill = celltype)) + 
  geom_violin(draw_quantiles = 0.5) + theme_bw() +
  geom_jitter(alpha = 0.5, width = 0.2, size = 3) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x = element_blank()) +
  ylab("Mean (logcounts)")
dev.off()

# scatterplot
plot.fname <- paste0("ggpt_",plot.fname.stem,".pdf")
pdf(file.path(plot.dpath, plot.fname), width = 5, height = 5)
ggplot(lct, aes(x = Neuron, y = `Non-neuron`, label = marker)) + 
  geom_point(alpha = 0.5, size = 3) + 
  geom_label_repel(box.padding   = 0.35, point.padding = 0.8, 
                   segment.color = 'grey50') + 
  theme_bw() + geom_abline(intercept = 0, slope = 1, col = "red")
dev.off()

