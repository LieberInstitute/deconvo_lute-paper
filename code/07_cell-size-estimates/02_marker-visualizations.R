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
# heatmap of marker logcounts
Heatmap(assays(setf)$logcounts)

# tile plot
lct <- assays(setf)$logcounts
lct$marker <- rownames(lct)
dfp <- rbind(data.frame(mean = lct$Neuron, marker = lct$marker),
             data.frame(mean = lct$`Non-neuron`, marker = lct$marker))
dfp$celltype <- c(rep("Neuron", nrow(lct)), rep("Non-neuron", nrow(lct)))
dfp$mean <- round(dfp$mean, 2)
ggplot(dfp, aes(x = celltype, y = marker, label = mean, fill = mean)) + 
  geom_tile() + geom_text() + theme_bw()

# violin plot
ggplot(dfp, aes(x = celltype, y = mean, fill = celltype)) + 
  geom_violin(draw_quantiles = 0.5) + theme_bw()

# violin with scatter
ggplot(dfp, aes(x = celltype, y = mean, fill = celltype)) + 
  geom_violin(draw_quantiles = 0.5) + theme_bw() +
  geom_jitter(alpha = 0.5, width = 0.2, size = 5)

# scatterplot
ggplot(lct, aes(x = Neuron, y = `Non-neuron`, label = marker)) + 
  geom_point(alpha = 0.5, size = 5) + 
  geom_label_repel(box.padding   = 0.35, point.padding = 0.8, 
                   segment.color = 'grey50') + 
  theme_bw() + geom_abline(intercept = 0, slope = 1, col = "red")
