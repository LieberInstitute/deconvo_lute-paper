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

trueProportionsVector <- seq(0.4, 0.9, 0.035)
bulkExpressionVector <- rnorm(
  length(trueProportionsVector), mean = 0.7, sd = 0.09)

listScatterResult <- multiPanelScatterPlots(
  trueProportionValueVector = trueProportionsVector,
  bulkExpressionVector = bulkExpressionVector
)

#--------------------------------------------------

# save

#--------------------------------------------------

#
#
#
#
# save multipanel scatterplots

jpeg(paste0("./figures/09_example_plots/",
            "scatterplot_samples_true-versus-predicted.jpg"), 
     width = 10, height = 2.5, units = "in", res = 400)

listScatterResult$ggScatter + xlab('Known')

dev.off()

#
#
#
#
# save image
save.image("./env/09_example_plots/03_scatterplot_multi_condition_script.RData")