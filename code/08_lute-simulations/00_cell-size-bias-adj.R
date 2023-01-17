#!/usr/bin/env R

# Author: Sean Maden
#
# Cell size bias adjustment experiments
#

libv <- "lute"
sapply(libv, library, character.only = TRUE)

#----------------
# run simulations
#----------------
lexpt <- donor_marker_biasexpt(offsetv = c(1, 10), P = c(0.25, 0.75),
                               donor.adj.method = 'var_denom',
                               gindexv = c(1, 2), ndonor = 10, ktotal = 2,
                               seed.num = 0, verbose = FALSE)