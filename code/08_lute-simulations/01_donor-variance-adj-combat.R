#!/usr/bin/env R

# Author: Sean Maden
#
# Show results of testing ComBat adjustment with varying donor-specific biases.
#

libv <- c("lute", "ggpubr", "gridExtra")
sapply(libv, library, character.only = TRUE)

# new save dest
save.dpath <- file.path("deconvo_method-paper", "outputs", "08_lute-simulations")

#----------------------
# get bias expt results
#----------------------
lb <- donor_marker_biasexpt(offsetv = c(5, 100), 
                            donor.adj.method = "combat", 
                            verbose = T)

lb$dfres