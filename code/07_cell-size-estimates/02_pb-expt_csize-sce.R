#!/usr/bin/env R

# Author: Sean Maden
#
# Show the impact of including cell size adjustments on pseudobulking outcomes.

libv <- c("lute", "SingleCellExperiment", "SummarizedExperiment", "ggplot2")
sapply(libv, library, character.only = T)

#----------
# load data 
#----------
save.dpath <- file.path("deconvo_method-paper", "outputs", "07_cell-size-estimates")

# cell size data
read.dpath <- file.path("deconvo_method-paper", "outputs", "07_cell-size-estimates")
sce.csize.fname <- "df-cellsize_donor-region_sce.rda"
df.csize <- get(load(file.path(read.dpath, sce.csize.fname)))

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
dfs <- aggregate(df.csize, by = list(df.csize$celltype), FUN = mean)
dfs <- dfs[,c(1,3)]; colnames(dfs) <- c("celltype", "mean_size")

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

#-------------------------------
# visualizations -- sim outcomes
#-------------------------------
plot.fnstem <- "lutesim-stransform_k2-n20-perk"
plot.dpath <- file.path("deconvo_method-paper", "outputs", "07_cell-size-estimates")

# facet scatter plot
plot.fname <- paste0("ggpt-facet-rmse_",plot.fnstem,".pdf")
pdf(file.path(plot.dpath, plot.fname), 5, 3)
lres$lgg$ggpt1 + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
dev.off()

# violin plot
plot.fname <- paste0("ggvp-rmse_",plot.fnstem,".pdf")
pdf(file.path(plot.dpath, plot.fname), 2, 2.5)
lres$lgg$ggvp + theme_bw() +
  geom_jitter(alpha = 0.1, width = 0.1, size = 0.1,
              col = "black") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Transformation")
dev.off()

# scatterplot
plot.fname <- paste0("ggpt-rmse_",plot.fnstem,".pdf")
pdf(file.path(plot.dpath, plot.fname), 5, 5)
lres$lgg$ggpt2 + geom_point(alpha = 0.5) + 
  theme_bw() + ggtitle("RMSE") + xlab("FALSE") + ylab("TRUE")
dev.off()