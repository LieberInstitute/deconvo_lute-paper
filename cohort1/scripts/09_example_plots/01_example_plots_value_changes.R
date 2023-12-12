#!/usr/bin/env R

libv <- c("ggplot2", "reshape2", "cowplot")
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


#
#
#
#
#
# save multipanel plot with shared legend, separated columns by variable type
plotLegend <- get_legend(listMultiPlot05$ggMulti)

dfPlotAll <- listMultiPlot05$dfPlotAll
dfPlotAll1 <- dfPlotAll[dfPlotAll$variableType=="condition",]
dfPlotAll2 <- dfPlotAll[dfPlotAll$variableType=="result",]


ggMultiPanel1 <- 
  ggplot(dfPlotAll1, aes(x = variable, y = value, fill = Change, color = variableType)) + 
  geom_bar(stat="identity", size = 1.2) + theme_bw() +
  ylab("Change (New - Old)") + facet_wrap(~conditionLabel, nrow = 1) + 
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none",
        axis.title.y = element_blank(), axis.title.x = element_blank()) + 
  scale_fill_manual(breaks = c("Increase", "Decrease"), 
                    values=c("dodgerblue", "gold")) +
  scale_color_manual(breaks = c("condition", "result"), 
                     values=c("#44AA99", "#332288")) +
  geom_text(aes(label = value), vjust = ifelse(dfPlotAll1$value >= 0, -1, 1)) +
  ylim(min(dfPlotAll1$value)-1.5, max(dfPlotAll1$value)+1.5)

ggMultiPanel2 <- 
  ggplot(dfPlotAll2, aes(x = variable, y = value, fill = Change, color = variableType)) + 
  geom_bar(stat="identity", size = 1.2) + theme_bw() +
  ylab("Change (New - Old)") + facet_wrap(~conditionLabel, nrow = 1) + 
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none",
        axis.title.y = element_blank(), axis.title.x = element_blank()) + 
  scale_fill_manual(breaks = c("Increase", "Decrease"), 
                    values=c("dodgerblue", "gold")) +
  scale_color_manual(breaks = c("condition", "result"), 
                     values=c("#44AA99", "#332288")) +
  geom_text(aes(label = value), vjust = ifelse(dfPlotAll2$value >= 0, -1, 1)) +
  ylim(min(dfPlotAll2$value)-0.2, max(dfPlotAll2$value)+0.2)

ggMultiPanelLegend <- 
  ggplot(dfPlotAll, 
         aes(x = variable, y = value, fill = Change)) + 
  geom_bar(stat="identity", size = 1.2) + theme_bw() +
  ylab("Change (New - Old)") + facet_wrap(~conditionLabel, nrow = 1) + 
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Affect of scale change") +
  scale_fill_manual(breaks = c("Increase", "Decrease"), 
                    values=c("dodgerblue", "gold")) +
  geom_text(aes(label = value), 
            vjust = ifelse(dfPlotAll$value >= 0, -1, 1)) +
  ylim(min(dfPlotAll$value)-1.5, max(dfPlotAll$value)+1.5)
ggMultiPanelLegend <- get_legend(ggMultiPanelLegend)

jpeg(paste0("./figures/09_example_plots/",
            "barplots_scale-changes_conditions-and-results-separate.jpg"), 
     width = 10.5, height = 6, units = "in", res = 400)
grid.arrange(ggMultiPanel1, ggMultiPanel2, ggMultiPanelLegend,
             layout_matrix = matrix(c(rep(1,5),3,rep(2,5),3), nrow = 2, byrow=T),
             left = paste0(paste0(rep(" ",10), collapse=""), "Value change (New - Old)"),
             bottom = paste0("Variable", paste0(rep(" ",30), collapse="")),
             top = paste0("Affect of scale change", paste0(rep(" ",170), collapse="")))
dev.off()

#
#
#
#
#
# save image
save.image(
  "./env/09_example_plots/01_example_plots_value_changes_script.RData")