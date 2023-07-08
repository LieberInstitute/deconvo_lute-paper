#!/usr/bin/env R

# Author: Sean Maden
#
# Perform constrained simulations with lute using the following params:
# K = 2 types total
# S: s1 > s2 (with slight random variation)
# P: p1~0.75; p2~0.25 (with slight random variation)
#
#

libv <- c("lute", "SingleCellExperiment", "SummarizedExperiment",
          "ComplexHeatmap", "ggplot2", "ggrepel")
sapply(libv, library, character.only = T)

#----------
# load data 
#----------
save.dpath <- file.path("deconvo_method-paper", "outputs", "07_cell-size-estimates")

# cell size data
csize.fname <- "dfcellsize_halo.rda"
csize.fpath <- file.path("deconvo_method-paper", "outputs", 
                         "07_cell-size-estimates", csize.fname)
dfs <- get(load(file.path(csize.fpath)))

# se marker data, k2
sef.fname <- "sef_mr-markers_k2_20-per-k_dlpfc-ro1.rda"
sef.dpath <- file.path("deconvo_method-paper", "outputs", 
                        "05_marker-gene-annotations")
sef <- get(load(file.path(sef.dpath, sef.fname)))

#-------------------------
# simulations -- sce.csize
#-------------------------
# get signature matrix, z
sef[["donor"]] <- sef[["BrNum"]]
setf <- set_from_sce(sef, groupvar = "donor", method = "mean",
                     assayname = "logcounts")
lct <- assays(setf)$logcounts

# get cell sizes, s
dfs$celltype <- ifelse(dfs$type %in% c("Excit", "Inhib"), "Neuron", "Non-neuron")
dfs <- aggregate(dfs, by = list("celltype"), FUN = mean)

# run simulations
num.sim <- 1e3
# make lgv
lgv <- lapply(seq(num.sim), function(ii){list(lct[,1], lct[,2])})
# make lpv
lpv <- make_lpv()
# make ls
dfs$k2 <- ifelse(grepl("Excit|Inhib", dfs$celltype), "Neuron", "Non-neuron")
lsv <- lapply(seq(num.sim), function(ii){
  c(mean(dfs[dfs$k2 == "Neuron",2]), mean(dfs[dfs$k2 == "Non-neuron",2]))
})

# run sims
lres <- decon_analysis(lgv = lgv, lpv = lpv, lsv = lsv)

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

#-------------------------------
# visualizations -- sim outcomes
#-------------------------------
# facet scatter plot
lres$lgg$ggpt1

# violin plot
lres$lgg$ggvp

# scatterplot
lres$lgg$ggpt2
