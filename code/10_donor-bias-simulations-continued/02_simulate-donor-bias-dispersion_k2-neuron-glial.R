#!/usr/bin/env R

# Author: Sean Maden
# 
# Testing the donor bias experiment
#

libv <- c("lute", "ggplot2")
sapply(libv, library, character.only = TRUE)

#-------------------------
# simulate multiple donors
#-------------------------
simulate_sce_donor_bias <- function(donor.coeff.vector = NULL, ...){
  require(ggplot2)
  require(SummarizedExperiment)
  require(SingleCellExperiment)
  lr <- list()
  if(is(donor.coeff.vector, "NULL")){
    donor.coeff.vector <- rnorm(10, mean = 200, sd = 30)
  }
  sce <- do.call(cbind, lapply(donor.coeff.vector, function(ii){
    scei <- get_sce_donor_bias(ii, ...)
    scei[["donor.coeff"]] <- ii; scei
    # as(scei, "SummarizedExperiment")
  })) 
  # sce <- as(sce, "SingleCellExperiment")
  # plot
  dfp <- do.call(rbind, lapply(donor.coeff.vector, function(ii){
    ctf <- counts(sce[,sce[["donor.coeff"]]==ii])
    dfi <- data.frame(mean = rowMeans(ctf), var = rowVars(ctf))
    dfi$donor.coeff <- ii; dfi
  }))
  dfp$donor.coeff <- as.character(dfp$donor.coeff)
  ggsm <- ggplot(dfp, aes(x = mean, y = var, color = donor.coeff)) + 
    geom_abline(slope = 1, intercept = 0) + 
    scale_x_log10() + scale_y_log10() + geom_smooth() +
    theme(legend.position = 'none') +
    xlab('Mean (log10)') + ylab("Var (log10)")
  # get the return object
  lr[["sce"]] <- sce; lr[["dfp"]] <- dfp; lr[["ggsmooth"]] <- ggsm
  return(lr)
}

ld <- simulate_sce_donor_bias()

ld$ggsmooth

ld$ggsmooth + facet_wrap(~donor.coeff)

#------------------------------------------------
# simulate donor bias for neurons and glial cells
#------------------------------------------------
# main params
set.seed(0)
num.genes <- 1000
num.neuron <- 1000
num.glial <- 1000
dc.sd.neuron <- 30
dc.mean.neuron <- 100
dc.sd.glial <- 80
dc.mean.glial <- 500

# neurons 
mean.vector <- seq(1, 40, 1e-3)
mean.vector <- sample(seq(length(mean.vector)), size = num.genes)
donor.coeff.vector <- rnorm(n = num.neuron, mean = dc.mean.neuron, sd = dc.sd.neuron)
ld.neuron <- simulate_sce_donor_bias(donor.coeff.vector = donor.coeff.vector,
                                     mean.vector = mean.vector)

# glial cells
mean.vector <- seq(1, 30, 1e-3)
mean.vector <- sample(seq(length(mean.vector)), size = num.genes)
donor.coeff.vector <- rnorm(n = num.glial, mean = dc.mean.glial, sd = dc.sd.glial)
ld.glial <- simulate_sce_donor_bias(donor.coeff.vector = donor.coeff.vector,
                                     mean.vector = mean.vector)

#--------------------------------
# setup : test donor bias effects
#--------------------------------
# expect donors to show differences in dispersion
# gindexv <- c(rep(1, 1000), rep(2, 1000))
# ddf <- random_donordf(ndonor = 10, gindexv = gindexv)
# dim(ddf) # [1] 4000   15

set.seed(1)
# set number of donors
ndonor <- 10
# set dispersions for set of random donors
gamma.pos <- rnorm(ndonor, mean = 10, sd = 5)
gamma.neg <- rnorm(ndonor, mean = 10, sd = 5)
gindexv <- c(rep(1, 1000), rep(2, 1000))
ddf <- do.call(rbind, lapply(seq(ndonor), function(ii){
  gp <- gamma.pos[ii]; gn <- gamma.neg[ii]
  ddf <- random_donordf(ndonor = 1, gindexv = gindexv, 
                        gamma.pos = gp, gamma.neg = gn)
  ddf$gamma.pos = gp; ddf$gamma.neg = gn
  ddf$donor.id <- paste0("donor", ii)
  return(ddf)
}))
colnames(ddf)[1] <- "counts"

# plot marker densities
cnv <- colnames(ddf)
unique.donors <- cnv[grepl("donor[0-9]*", cnv)]
dfp <- do.call(rbind, lapply(unique.donors, function(donori){
  dfi <- data.frame(counts = ddf[,donori], 
                    type = ddf$type, 
                    marker.type = ddf$marker.type)
  dfi$donor <- donori; return(dfi)
}))


ggdens <- ggplot(ddf, aes(x = counts, color = donor)) + 
  geom_density(position = "stack") + theme_bw()

ggdens

ggdens + facet_wrap(~type)

ggdens + facet_wrap(~marker.type)
