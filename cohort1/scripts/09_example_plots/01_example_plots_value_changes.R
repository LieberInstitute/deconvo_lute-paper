#!/usr/bin/env R

libv <- c("ggplot2", "reshape2")
sapply(libv, library, character.only = T)

# Author: Sean Maden
#
# Plot to show direction change in expression and prediction with scale factor change.
#
#
#
#
#

source("./scripts/09_example_plots/00_param.R")

#--------------------------------------------------

# multiple panels -- marker expression start is 0.5

#--------------------------------------------------
# get multiple panels
# view facet of multi panel plots

# marker expression start is 0.5
listMultiPlot05 <- multiPanelPlots(markerExpressionStart = 0.5)
barplotStart05 <- 
  listMultiPlot05$resultsList$Decrease$result$valuesList$ggBarplotStart
grid.arrange(barplotStart05, listMultiPlot05$ggMulti, nrow = 1,
             layout_matrix = matrix(c(1,2,2,2),nrow=1))
barplotStart05WithLegend <- 
  barplotStart05 + theme(legend.position = "right")

# marker expression start is 1
listMultiPlot1 <- multiPanelPlots(markerExpressionStart = 1)
barplotStart1 <- 
  listMultiPlot1$resultsList$Decrease$result$valuesList$ggBarplotStart
grid.arrange(barplotStart1, listMultiPlot1$ggMulti, nrow = 1,
             layout_matrix = matrix(c(1,2,2,2),nrow=1))

# marker expression start is 2
listMultiPlot2 <- multiPanelPlots(markerExpressionStart = 2)
barplotStart2 <- 
  listMultiPlot2$resultsList$Decrease$result$valuesList$ggBarplotStart
grid.arrange(barplotStart2, listMultiPlot2$ggMulti, nrow = 1,
             layout_matrix = matrix(c(1,2,2,2),nrow=1))

#--------------------------

# save

#--------------------------

# save single panel plots
# note: save this once, because it does not change across conditions
jpeg(paste0("./figures/09_example_plots/",
            "scalefactors_example_barplots_marker05.jpg"), width = 5, 
     height = 3, units = "in", res = 400)
listMultiPlot05$ggBarCellScaleFactors
dev.off()

#
#
#
#
#
# save point values for simulation start barplot
# without legend
ggsave(filename = 
         paste0("./figures/09_example_plots/",
                "sim_start_expr05_barplot.jpg"),
       plot = barplotStart05,
       device = "jpeg", width = 3, height = 3, 
       units = "in", dpi = 400)
# with legend
ggsave(filename = 
         paste0("./figures/09_example_plots/",
                "sim_start_expr05_barplot-with-legend.jpg"),
       plot = barplotStart05WithLegend,
       device = "jpeg", width = 3.5, height = 3, 
       units = "in", dpi = 400)

# save multipanel change plots

# starting marker expression 05
ggsave(filename = 
         paste0("./figures/09_example_plots/",
                "multipanel_change_expr05_barplots.jpg"),
       plot = listMultiPlot05$ggMulti,
       device = "jpeg", width = 12, height = 4, 
       units = "in", dpi = 400)

# starting marker expression 1
ggsave(filename = 
         paste0("./figures/09_example_plots/",
                "multipanel_change_expr1_barplots.jpg"),
       plot = listMultiPlot05$ggMulti,
       device = "jpeg", width = 12, height = 4, 
       units = "in", dpi = 400)


# save marker expression barplots

# save multipanel plots
# save
jpeg(paste0("./figures/09_example_plots/",
            "multipanel_example_barplots_marker05.jpg"), width = 11, 
     height = 4, units = "in", res = 400)
grid.arrange(barplotStart05, listMultiPlot05$ggMulti, nrow = 1,
             layout_matrix = matrix(c(1,2,2,2),nrow=1))
dev.off()

jpeg(paste0("./figures/09_example_plots/",
            "multipanel_example_barplots_marker1.jpg"), width = 11, 
     height = 4, units = "in", res = 400)
grid.arrange(barplotStart1, listMultiPlot1$ggMulti, nrow = 1,
             layout_matrix = matrix(c(1,2,2,2),nrow=1))
dev.off()

jpeg(paste0("./figures/09_example_plots/",
            "multipanel_example_barplots_marker2.jpg"), width = 11, 
     height = 4, units = "in", res = 400)
grid.arrange(barplotStart2, listMultiPlot2$ggMulti, nrow = 1,
             layout_matrix = matrix(c(1,2,2,2),nrow=1))
dev.off()

# save image
save.image(
  "./env/09_example_plots/01_example_plots_value_changes_script.RData")