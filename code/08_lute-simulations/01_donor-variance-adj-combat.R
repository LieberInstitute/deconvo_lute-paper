#!/usr/bin/env R

# Author: Sean Maden
#
# Show results of testing ComBat adjustment with varying donor-specific biases.
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