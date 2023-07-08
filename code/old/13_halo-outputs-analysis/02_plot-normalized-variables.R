#!/usr/bin/env R

# Author: Sean Maden
#
# Plot summaries of normalized variables for cell area from HALO outputs. Note, 
# log10 transformation doesn't increase normality in AKT3 marker data, so we 
# leave it out of ANOVA tests.
#

source("deconvo_method-paper/code/_halo-outputs-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
halo.outputs.table <- get(load(output.updated.path))
halo.outputs.table <- halo.outputs.table %>% as.data.frame()

# plot nucleus area
# prepare plot data
plot.data1 <- data.frame(value = halo.output.table[,cell.area.variable])
plot.data2 <- data.frame(value = halo.output.table[,normalized.area.variable])
plot.data1$variable <- "original"
plot.data2$variable <- "log10"
plot.data.final <- rbind(plot.data1, plot.data2)
# save
jpeg(boxplot.area.jpg.path, width = 3.5, height = 3, units = "in", res = 400)
ggplot(plot.data.final, aes(x = variable, y = value)) + geom_boxplot() + 
  facet_zoom(ylim = c(0, 6)) + xlab("Normalization") + ggtitle("Nucleus area") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# plot marker
# prepare plot data
plot.data1 <- data.frame(value = halo.output.table[,gene.marker.label])
plot.data2 <- data.frame(value = halo.output.table[,normalized.marker.variable])
plot.data1$variable <- "original"
plot.data2$variable <- "log10"
plot.data.final <- rbind(plot.data1, plot.data2)
# save
jpeg(boxplot.marker.jpg.path, width = 3.5, height = 3, units = "in", res = 400)
ggplot(plot.data.final, aes(x = variable, y = value)) + geom_boxplot() + 
  facet_zoom(ylim = c(0, 6)) + xlab("Normalization") + ggtitle("AKT3 expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()