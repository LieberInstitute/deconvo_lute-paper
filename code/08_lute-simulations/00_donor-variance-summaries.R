#!/usr/bin/env R

# Author: Sean Maden
#
# Summarize results of varying the donor-specific variances.
#

libv <- c("lute", "ggpubr", "gridExtra")
sapply(libv, library, character.only = TRUE)

# new save dest
save.dpath <- file.path("deconvo_method-paper", "outputs", "08_lute-simulations")

#------------------------------------------
# simulate low and high variance donor data
#------------------------------------------
# show impact of increasing donor offset
ndonor <- 10
gindexv <- c(1, 2)
offset.low <- 1
offset.high <- 100

# get pca results
ldat <- list()

# simulate donor data
# low variances
title.str <- paste0("Low donor variance (offset = ",1,")\n")
ldat[["low"]] <- rand_donor_marker_table(ndonor = ndonor, gindexv = gindexv, 
                               sd.offset.pos = offset.low, 
                               sd.offset.neg = offset.low)
# high variances
ldat[["high"]] <- rand_donor_marker_table(ndonor = ndonor, gindexv = gindexv, 
                               sd.offset.pos = offset.high, 
                               sd.offset.neg = offset.high)

#------------------------------
# get PCA results of donor expt
#------------------------------
lpca <- list()
lpca[["low"]] <- pcaplots_donor(dt = ldat[["low"]], title.append = title.str)
lpca[["high"]] <- pcaplots_donor(dt = ldat[["high"]], title.append = title.str)

# save new results object
ldat[["pca_results"]] <- lpca
save.fname <- "lpca-results_donorvar-low-high_lute-donorsim.rda"
save(ldat, file = file.path(save.dpath, save.fname))

#-----------------
# plot pca results
#-----------------
# get scatterplot objects
# get scatterplot plot legend
ggpt1.leg <- lpca$low$pca.bydonortype$scatterplot.pc1.pc2 +
  scale_shape_discrete(labels = c("Neuron", "Non-neuron"),
                       name = "celltype") +
  ggtitle("Low donor variance")
ggleg <- get_legend(ggpt1.leg)
# get final scatterplots
ggpt2.leg <- lpca$high$pca.bydonortype$scatterplot.pc1.pc2 +
  scale_shape_discrete(labels = c("Neuron", "Non-neuron"),
                       name = "celltype") +
  ggtitle("High donor variance")
ggpt2.noleg <- ggpt2.leg + theme(legend.position = 'none')
ggpt1.noleg <- ggpt1.leg + theme(legend.position = 'none')

# save pca scatterplots
# save composite
lm <- matrix(c(1,1,2,2,3), nrow = 1)
plot.fname <- "pca-composite_donorvar-low-high_lute-sim.jpg"
jpeg(file.path(save.dpath, plot.fname), width = 9.7, height = 4.2, 
     units = "in", res = 400)
grid.arrange(ggpt1.noleg, ggpt2.noleg, ggleg, layout_matrix = lm)
dev.off()
# low var
plot.fname <- "pca_donorvar-low_lute-sim.jpg"
jpeg(file.path(save.dpath, plot.fname), width = 7, height = 5, 
    units = "in", res = 400)
ggpt1.leg
dev.off()
# high var
plot.fname <- "pca_donorvar-high_lute-sim.jpg"
jpeg(file.path(save.dpath, plot.fname), width = 7, height = 5, 
    units = "in", res = 400)
ggpt2.leg
dev.off()
