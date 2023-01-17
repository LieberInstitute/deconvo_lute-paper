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
# define offset series
offsetv <- c(seq(0.1, 2.1, 0.1), seq(18, 20, 0.1))


lexpt <- donor_marker_biasexpt(offsetv = offsetv, P = c(0.25, 0.75),
                               donor.adj.method = 'var_denom',
                               gindexv = c(1, 2), ndonor = 10,
                               seed.num = 0, verbose = FALSE)