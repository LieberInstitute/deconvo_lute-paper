#!/usr/bin/env R

# Author: Sean Maden
# 
# Validate sopt utility on new held-out bulk samples.
#

libv <- c("lute")
sapply(libv, library, character.only = T)

# load
# validation data