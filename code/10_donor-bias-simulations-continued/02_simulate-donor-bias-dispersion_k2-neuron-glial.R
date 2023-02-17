#!/usr/bin/env R

# Author: Sean Maden
# 
# Testing the donor bias experiment
#

libv <- c("lute", "ggplot2")
sapply(libv, library, character.only = TRUE)


#------------------
# get splatter data
#------------------
BiocManager::install("splatter")
library(splatter)

params <- splatEstimate(sce)
sim <- splatSimulate(params)


# Load package
suppressPackageStartupMessages({
  library(splatter)
  library(scater)
})

# Create mock data
set.seed(1)
sce <- mockSCE()

# Estimate parameters from mock data
params <- splatEstimate(sce)
#> NOTE: Library sizes have been found to be normally distributed instead of log-normal. You may want to check this is correct.
# Simulate data using estimated parameters
sim <- splatSimulate(params)

plot_mean_var_sce <- function(sce, donor.variable = NULL){
  if(is(donor.variable, "NULL")){
    ct <- counts(sce)
    dfp <- data.frame(mean = rowMeans(ct), var = rowVars(ct))
    ggplot(dfp, aes(x = mean, y = var)) +
      geom_smooth() + geom_abline(intercept = 0, slope = 1) +
      theme_bw() + scale_x_log10() + scale_y_log10()
  } else{
    donor.vector <- unique(sce[[donor.variable]])
    dfp <- do.call(rbind, lapply(donor.vector, function(donori){
      ctf <- counts(sce[,sce[[donor.variable]]==donori])
      dfi <- data.frame(mean = rowMeans(ctf), var = rowVars(ctf))
      dfi$donor <- donori; dfi
    }))
    dfp$donor <- as.character(dfp$donor)
    ggplot(dfp, aes(x = mean, y = var, color = donor)) +
      geom_smooth() + geom_abline(intercept = 0, slope = 1) +
      theme_bw() + scale_x_log10() + scale_y_log10()
  }
}

plot_mean_var_sce(sce)
plot_mean_var_sce(sim)

donors_from_sce <- function(sce, bias.vector = seq(3), bias.means = seq(3)){
  ct <- counts(sce)
  medians <- rowMedians(ct) + 1
  num.cells <- ncol(ct)
  num.genes <- nrow(ct)
  sce.new <- do.call(cbind, lapply(bias.vector, function(biasi){
    sce.new <- sce
    d <- biasi*(medians/(medians^2))
    #biasv <- unlist(sapply(d, function(di){
    #  rnbinom(1, size = di, mu = 1)
    #}))
    #ct.new <- sweep(ct, 1, biasv, FUN = "+")
    #counts(sce.new) <- ct.new
    ct.new <- do.call(rbind, lapply(seq(num.genes), function(ii){
      rnbinom(num.cells, size = d[ii], mu = medians[ii])
    }))
    rownames(ct.new) <- rownames(ct)
    colnames(ct.new) <- colnames(ct)
    assays(sce.new)[["counts"]] <- ct.new
    sce.new[["donor.bias"]] <- biasi
    return(sce.new)
  }))
  return(sce.new)
}

set.seed(0)
sce1 <- mockSCE()
plot_mean_var_sce(sce1)

params <- splatEstimate(sce1)
sce2 <- splatSimulate(params)
plot_mean_var_sce(sce2)

sce3 <- donors_from_sce(sce2, c(100, 200, 300), rep(5, 3))
plot_mean_var_sce(sce3, "donor.bias")

plot_mean_var_sce(sce3, "donor.bias") + facet_wrap(~donor)

#-------------------------
# simulate multiple donors
#-------------------------

get_sce_donor_bias <- function(donor.bias.coeff = 1, mean.vector = NULL,
                               total.cells = 10, num.genes.iter = 10,
                               num.types = 1, seed.num = 0){
  set.seed(seed.num)
  if(is(mean.vector, "NULL")){mean.vector = seq(1, 30, 1)}
  sce <- do.call(rbind, lapply(mean.vector, function(meani){
    d <- donor.bias.coeff*(meani/(meani^2))
    scei <- random_sce(num.types = num.types, num.cells = total.cells, 
                       expr.mean = meani, num.genes = num.genes.iter, 
                       dispersion = d)
    # as(scei, "SummarizedExperiment")
    return(scei)
  }))
  return(as(sce, "SingleCellExperiment"))
}

simulate_sce_donor_bias <- function(donor.coeff.vector = NULL,
                                    mean.vector = NULL, 
                                    cells.per.donor = 1000,
                                    num.genes.iter = 10){
  lr <- list()
  if(is(donor.coeff.vector, "NULL")){
    donor.coeff.vector <- rnorm(10, mean = 200, sd = 30)
  }
  sce <- do.call(cbind, lapply(donor.coeff.vector, function(ii){
    scei <- get_sce_donor_bias(donor.bias.coeff = ii, 
                               mean.vector = mean.vector, 
                               total.cells = cells.per.donor,
                               num.genes.iter = num.genes.iter)
    scei[["donor.coeff"]] <- ii; scei
  }))
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
# mean vector
num.mean.bin <- 100
# iteration params
num.genes.iter <- 10
cells.per.donor <- 500
# donor coeff params
num.donor <- 10

# neurons 
dc.sd.neuron <- 25
dc.mean.neuron <- 80
max.mean.neuron <- 100
# mean.vector <- seq(1, max.mean.neuron, 1e-3)
mean.sample <- sample(seq(length(mean.vector)), size = num.mean.bin)
#mean.vector <- mean.vector[mean.sample]
mean.vector <- rnbinom(n = 10000, size = 1e-1, mu = 0.5)
mean.vector <- mean.vector[mean.sample]
plot(density(mean.vector))
summary(mean.vector)
mean.vector.neuron <- mean.vector
donor.coeff.vector <- rnorm(n = num.donor, mean = dc.mean.neuron, 
                            sd = dc.sd.neuron)
ld.neuron <- simulate_sce_donor_bias(donor.coeff.vector = donor.coeff.vector,
                                     mean.vector = mean.vector, 
                                     cells.per.donor = cells.per.donor,
                                     num.genes.iter = num.genes.iter)
ld.neuron[["mean.vector"]] <- mean.vector.neuron

# glial cells
max.mean.glial <- 35
dc.sd.glial <- 80
dc.mean.glial <- 300
# mean.vector <- seq(1, max.mean.glial, 1e-3)
mean.sample <- sample(seq(length(mean.vector)), size = num.mean.bin)
# mean.vector <- mean.vector[mean.sample]
mean.vector <- rnbinom(n = 10000, size = 9e-3, mu = 1.1e-1)
mean.vector <- mean.vector[mean.sample]
plot(density(mean.vector))
summary(mean.vector)
mean.vector.glial <- mean.vector

donor.coeff.vector <- rnorm(n = num.donor, mean = dc.mean.glial, sd = dc.sd.glial)
ld.glial <- simulate_sce_donor_bias(donor.coeff.vector = donor.coeff.vector,
                                    mean.vector = mean.vector,
                                    cells.per.donor = cells.per.donor,
                                    num.genes.iter = num.genes.iter)
ld.glial[["mean.vector"]] <- mean.vector.glial

# save results
ld <- list(neuron = ld.neuron, glial = ld.glial)
fname <- "lute-sim-donorbias_k2-sce-sim_ro1-dlpfc.rda"
save(ld, file = fname)

#---------------------
# plot simulation data
#---------------------
# plots
# plot neurons
ld.neuron$ggsmooth
ld.neuron$ggsmooth + facet_wrap(~donor.coeff)
# plot glial
ld.glial$ggsmooth
ld.glial$ggsmooth + facet_wrap(~donor.coeff)

# save composite plot
dfp1 <- ld.neuron$dfp
dfp1$type <- "neuron"
dfp2 <- ld.glial$dfp
dfp2$type <- "glial"
dfp <- rbind(dfp1, dfp2)
ggsm <- ggplot(dfp, aes(x = mean, y = var, color = donor.coeff)) + 
  geom_abline(slope = 1, intercept = 0) + theme_bw() +
  scale_x_log10() + scale_y_log10() + geom_smooth() +
  theme(legend.position = 'none') +
  xlab('Mean (log10)') + ylab("Var (log10)")
# save jpg
fname <- "ggsm-composite_lute-sim-donorbias_k2-mean-var_ro1-dlpfc.jpg"
jpeg(file = fname, width = 6, height = 3, units = "in", res = 400)
ggsm + facet_wrap(~type); dev.off()



#-------------------------------------------
# get composite with empirical/observed data
#-------------------------------------------
markeri <- "k2"
handle.str <- "ro1-dlpfc"
# load
# get save dpath
code.dname <- "10_donor-bias-simulations-continued"
proj.dname <- "deconvo_method-paper"
save.dpath <- file.path(proj.dname, "outputs", code.dname)
# variances
fname <- paste0("matrix-donor-type_variances_",markeri, "_", handle.str,".rda")
mvar <- get(load(file.path(save.dpath, fname)))
# means
fname <- paste0("matrix-donor-type_means_",markeri, "_", handle.str,".rda")
mmean <- get(load(file.path(save.dpath, fname)))

# get plot data
genev <- rownames(sce.mean)
groupv <- colnames(sce.mean)
dfp.obs <- do.call(rbind, lapply(seq(ncol(sce.mean)), function(ii){
  dfpi <- data.frame(gene = genev, mean = sce.mean[,ii], var = sce.var[,ii])
  dfpi$group <- groupv[ii]; return(dfpi)
}))
dfp.obs$celltype <- gsub(";.*", "", dfp.obs$group)
dfp.obs$donor <- gsub(".*;", "", dfp.obs$group)
dfp.obs$dataset <- "observed"


# bind sim data
# save composite plot
dfp1 <- ld.neuron$dfp
dfp1$type <- "neuron"
dfp2 <- ld.glial$dfp
dfp2$type <- "glial"
dfp.sim <- rbind(dfp1, dfp2)
dfp.sim$dataset <- "simulated"

# bind all
dfp.obs.filt <- dfp.obs[,c(2,3,6,5,7)]
colnames(dfp.sim) <- colnames(dfp.obs.filt)
dfp <- rbind(dfp.obs.filt, dfp.sim)

# plot
ggsm <- ggplot(dfp, aes(x = mean, y = var, color = donor)) + 
  geom_abline(slope = 1, intercept = 0) + theme_bw() +
  scale_x_log10() + scale_y_log10() + geom_smooth() +
  theme(legend.position = 'none') +
  xlab('Mean (log10)') + ylab("Var (log10)")

ggsm + facet_wrap(~dataset+celltype)

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
