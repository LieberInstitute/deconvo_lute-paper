#!/usr/bin/env R

# Author: Sean Maden
#
# Test independent pseudobulk from multi-region brain data.

source("deconvo_method-paper/code/14_deconvolution-framework-trials/00_parameters.R")
sapply(libv, library, character.only = T)
sce.mrb <- get(load(sce.mrb.path))