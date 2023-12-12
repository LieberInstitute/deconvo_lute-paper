#!/usr/bin/env R

libv <- c("ggplot2", "reshape2", "gridExtra", "cowplot")
sapply(libv, library, character.only = T)

# Author: Sean Maden
#
# Gets scatterplots of multiple conditions, varying the true proportions and bulk expressions.
#
#
#
#

#------
# load
#------
source("./scripts/09_example_plots/00_param.R")

load("./env/09_example_plots/01_example_plots_value_changes_script.RData")

#--------------------------------------------------

# multiple panels -- marker expression start is 0.5

#--------------------------------------------------
# get multiple panels
# view facet of multi panel plots

# marker expression start is 0.5
listMultiPlot05 <- multiPanelPlots(markerExpressionStart = 0.5)


