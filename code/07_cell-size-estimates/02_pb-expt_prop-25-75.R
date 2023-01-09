#!/usr/bin/env R

# Author: Sean Maden
#
#

libv <- c("lute", "SingleCellExperiment", "SummarizedExperiment", "ggplot2")
sapply(libv, library, character.only = T)

#----------
# load data 
#----------
save.dpath <- file.path("deconvo_method-paper", "outputs", 
                        "07_cell-size-estimates")

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
df.csize$celltype <- ifelse(grepl("Excit|Inhib", df.csize$celltype), 
                            "Neuron", "Non-neuron")
dfs <- aggregate(df.csize, by = list(df.csize$celltype), FUN = mean)
dfs <- dfs[,c(1,3)]; colnames(dfs) <- c("celltype", "mean_size")
dfs$k2 <- dfs[,1]
# make ls
lsv <- lapply(seq(num.sim), function(ii){
  c(mean(dfs[dfs$k2 == "Neuron",2]), mean(dfs[dfs$k2 == "Non-neuron",2]))
})

# simulations params
num.sim <- 10
# make lgv
lgv <- lapply(seq(num.sim), function(ii){list(lct[,1], lct[,2])})
# make lpv
p1v <- seq(0.70, 0.80, 0.1/num.sim); p2v <- 1-p1v
lpv <- lapply(seq(num.sim), function(ii){c(p1v[ii], p2v[ii])})

# run sims
lres <- decon_analysis(lgv = lgv, lpv = lpv, lsv = lsv)

#-----------------
# make scatterplot
#-----------------
