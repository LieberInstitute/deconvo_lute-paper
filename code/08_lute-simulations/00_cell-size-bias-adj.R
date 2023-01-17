#!/usr/bin/env R

# Author: Sean Maden
#
# Cell size bias adjustment experiments
#

libv <- c("lute", "ggpubr", "gridExtra")
sapply(libv, library, character.only = TRUE)

# new save dest
save.dpath <- file.path("deconvo_method-paper", "outputs", "08_lute-simulations")

#------------------------------
# pca of low and high donor var
#------------------------------
# show impact of increasing donor offset
ndonor <- 10
gindexv <- c(1, 2)
offset.low <- 1
offset.high <- 100

# get pca results
lpca <- list()
# pca, low
title.str <- paste0("Low donor variance (offset = ",1,")\n")
dfp <- rand_donor_marker_table(ndonor = ndonor, gindexv = gindexv, 
                               sd.offset.pos = offset.low, 
                               sd.offset.neg = offset.low)
lpca[["low"]] <- pcaplots_donor(dt = dfp, title.append = title.str)
# pca.high
title.str <- paste0("High donor variance (offset = ",10,")\n")
dfp <- rand_donor_marker_table(ndonor = ndonor, gindexv = gindexv, 
                               sd.offset.pos = offset.high, 
                               sd.offset.neg = offset.high)
lpca[["high"]] <- pcaplots_donor(dt = dfp, title.append = title.str)

# get scatterplot objects
# get scatterplot plot legend
ggpt <- lpca$low$pca.bydonortype$scatterplot.pc1.pc2
ggleg <- get_legend(ggpt)
# get final scatterplots
ggpt1 <- ggpt
ggpt2 <- lpca$high$pca.bydonortype$scatterplot.pc1.pc2

# save pca scatterplots
# save composite
lm <- matrix(c(1,1,2,2,3), nrow = 1)
plot.fname <- "pca_donorvar-low_lute-sim.jpg"
jpeg(file.path(save.dpath, plot.fname), width = 7, height = 5, 
     units = "in", res = 400)
grid.arrange(ggpt1, ggpt2, ggleg, layout_matrix = lm)
dev.off()
# low var
plot.fname <- "pca_donorvar-low_lute-sim.jpg"
jpeg(file.path(save.dpath, plot.fname), width = 7, height = 5, 
    units = "in", res = 400)
lpca$low$pca.bydonortype$scatterplot.pc1.pc2
dev.off()
# high var
plot.fname <- "pca_donorvar-high_lute-sim.jpg"
jpeg(file.path(save.dpath, plot.fname), width = 7, height = 5, 
    units = "in", res = 400)
lpca$high$pca.bydonortype$scatterplot.pc1.pc2
dev.off()

# save pca screeplots



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
