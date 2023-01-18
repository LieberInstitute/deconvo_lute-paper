#!/usr/bin/env R

# Author: Sean Maden
#
# Cell size bias adjustment experiments
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

#--------------------
# test using sd_denom
#--------------------
# define offset series
set.seed(0)
pdiffv <- seq(-2e-2, 2e-2, 2e-3)
ddiffv <- rnorm(n = 21, mean = 1e-2, sd = 1e-3)

lseries <- lapply(seq(pdiffv), function(ii){
  p1true <- 0.75 + pdiffv[ii]; p2true <- 1-p1true
  donor_marker_biasexpt(offsetv = offsetv[ii], 
                        P = c(p1true, p2true),
                        donor.adj.method = "sd_denom",
                        bounds_thresh = 50)
})

dfr <- do.call(rbind, lapply(lseries, function(ii){ii$dfres}))

# scatterplots of results
dfp <- dfr
dfp$celltype <- ifelse(dfp$type.index==1, "neuron", "non-neuron")
dfp$offset.type <- ifelse(dfp$offset < 10 , "low", "high")

ggpt <- ggplot(dfp, aes(x = prop.true, y = prop.pred, 
                color = celltype, shape = celltype)) + theme_bw() + 
  geom_point(alpha = 0.5) + 
  geom_abline(slope = 1, intercept = 0, color = "black") +
  theme()

ggpt + facet_wrap(~prop.type*offset.type)

#-----------------------------------------
# test using limma::removeBatchEffects()
#-----------------------------------------
set.seed(0)
pdiffv <- seq(-2e-2, 2e-2, 2e-3)
ddiffv <- rnorm(n = length(pdiffv), mean = 1e-2, sd = 1e-3)

lseries <- lapply(seq(pdiffv), function(ii){
  p1true <- 0.75 + pdiffv[ii]; p2true <- 1-p1true
  donor_marker_biasexpt(offsetv = offsetv[ii], P = c(p1true, p2true),
                        donor.adj.method = "limma", bounds_thresh = 50)
})

dfr <- do.call(rbind, lapply(lseries, function(ii){ii$dfres}))

# scatterplots of results
dfp <- dfr
dfp$celltype <- ifelse(dfp$type.index==1, "neuron", "non-neuron")
dfp$offset.type <- ifelse(dfp$offset < 10 , "low", "high")

ggpt <- ggplot(dfp, aes(x = prop.true, y = prop.pred, 
                        color = celltype, shape = celltype)) + theme_bw() + 
  geom_point(alpha = 0.5) + 
  geom_abline(slope = 1, intercept = 0, color = "black") +
  theme()

ggpt + facet_wrap(~prop.type*offset.type)

#-----------------------------------------
# test using sctransform::correct_counts()
#-----------------------------------------
# install.packages("sctransform")
