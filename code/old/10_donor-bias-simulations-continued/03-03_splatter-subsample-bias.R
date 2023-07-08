#!/usr/bin/env R

# Author: Sean Maden
#
# Uses splatter to evaluate sources of within-sample/batch bias effects.
#
#

libv <- c("lute")
sapply(libv, library, character.only = T)

#-------------------
# do lute simulation
#-------------------
