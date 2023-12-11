#!/usr/bin/env R

libv <- c("ggplot2", "reshape2", "gridExtra", "cowplot")
sapply(libv, library, character.only = T)

# Author: Sean Maden
#
# Plots the simulation values barplots
#
#
#
#

#------
# load
#------
source("./scripts/09_example_plots/00_param.R")
load("./env/09_example_plots/01_example_plots_value_changes_script.RData")

#-------------------
# get barplot values
#-------------------
# get iterations labels
dfIter <- listMultiPlot05$dfPlotAll
changeLabels <- unique(gsub("\\..*", "", rownames(dfIter)))
variablesVector <- c("cellScaleFactor", "markerExpressionScaled",
                  "predictedProportion", "bias", "error")

# get plots data list
lgg <- lapply(changeLabels, function(changeLabel){
  valuesPlot <- 
    listMultiPlot05$resultsList[[changeLabel]]$result$valuesList
  deconvoResultNew <- 
    listMultiPlot05$resultsList[[changeLabel]]$result$deconvoResult
  deconvoResultNew <- deconvoResultNew@predictionsTable
  cellScaleFactor <- valuesPlot$cellScaleFactorNew
  markerExpressionScaled <- valuesPlot$zrefNewScaled[1,1]
  trueProportion <- valuesPlot$trueProportionValue
  predictedProportion <- deconvoResultNew[1,1]
  biasValue <- trueProportion-predictedProportion
  errorValue <- abs(biasValue)
  dfPlot <- data.frame(
    values = c(cellScaleFactor, markerExpressionScaled, 
               predictedProportion, biasValue, errorValue),
    variables = variablesVector,
    variableType = c(rep("condition", 2), rep("result", 3))
  )
  dfPlot$variables <- 
    factor(dfPlot$variables, levels = variablesVector)
  dfPlot$type <- changeLabel
  dfPlot$value <- round(dfPlot$value, 2)
  newBarPlot <- 
    ggplot(dfPlot, aes(x = variables, y = value, color = variableType)) + 
    geom_abline(intercept = 0, slope = 0) +
    geom_text(aes(label = value), 
              vjust = ifelse(dfPlot$value >= 0, -1.2, 1.2)) +
    geom_bar(stat = "identity", size = 1.2) + theme_bw() +
    scale_color_manual(breaks = c("condition", "result"), 
                       values=c("#44AA99", "#332288")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(changeLabel) +
    ylim(min(dfPlot$value)-1.2, max(dfPlot$value)+1.2)
  
  return(list(dfPlot = dfPlot,
              plot = newBarPlot))
})

# bind together all plot dfs
dfPlotAll <- do.call(rbind, lapply(lgg, function(item){
  return(item$dfPlot)
}))



#--------------------------
#
#
# get single composite plot
#
#
#--------------------------


# Example, one legend, six plots

plot1 <- barplotStart05
plot2 <- lgg$plot[[1]] + theme(legend.position = "none")
plot3 <- lgg$plot[[2]] + theme(legend.position = "none")
plot4 <- lgg$plot[[3]] + theme(legend.position = "none")
plot5 <- lgg$plot[[4]] + theme(legend.position = "none")
plot6 <- lgg$plot[[5]] + theme(legend.position = "none")
plotLegend <- get_legend(lgg[[1]])

grid.arrange(plot1, plot2, plot3, plot4,
             plot5, plot6, plotLegend, 
             layout_matrix = matrix(c(1,2,3,4,5,6,7,7), nrow = 2))

# facet barplots
newFacetBarplots <- 
  ggplot(dfPlotAll, aes(x = variables, y = value, color = variableType)) + 
  geom_abline(intercept = 0, slope = 0) +
  geom_text(aes(label = value), 
            vjust = ifelse(dfPlotAll$value >= 0, -1.2, 1.2)) +
  geom_bar(stat = "identity", size = 1.2) + theme_bw() +
  scale_color_manual(breaks = c("condition", "result"), 
                     values=c("#44AA99", "#332288")) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1)) +
  ylim(min(dfPlotAll$values)-1.2, max(dfPlotAll$values)+1.2) +
  facet_wrap(~type)

ggsave(filename = 
         paste0("./figures/09_example_plots/",
                "multipanel_simulation_values_barplots.jpg"),
       plot = newFacetBarplots,
       device = "jpeg", width = 10, height = 5, units = "in", dpi = 400)



#--------
# save
#--------

save.image("./env/09_example_plots/02_simulation_values_barplots_script.RData")