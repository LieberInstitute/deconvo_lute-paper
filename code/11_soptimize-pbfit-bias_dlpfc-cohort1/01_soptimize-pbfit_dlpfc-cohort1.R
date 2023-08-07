#!/usr/bin/env R

#
# Get pbfit across k2 experiments
#

source("deconvo_method-paper/code/11_soptimize-pbfit-bias_dlpfc-cohort1/00_parameters-pseudobulk.R")
sapply(libv, library, character.only = T)