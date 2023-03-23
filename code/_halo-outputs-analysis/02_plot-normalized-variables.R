#!/usr/bin/env R

# Author: Sean Maden
#
# Parameters for HALO image analyses.
#

source("deconvo_method-paper/code/_halo-outputs-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
halo.outputs.table <- get(load(output.updated.path))
halo.outputs.table <- halo.outputs.table %>% as.data.frame()

# format data for plot
marker.vector.original <- halo.output.table[,gene.marker.label]
marker.vector.normalized1 <- halo.output.table[,normalization.variable1]
marker.vector.normalized2 <- halo.output.table[,normalization.variable2]
plot.data1 <- data.frame(value = marker.vector.original)
plot.data2 <- data.frame(value = marker.vector.normalized1)
plot.data3 <- data.frame(value = marker.vector.normalized2)
plot.data1$variable <- "original"
plot.data2$variable <- "normalization1"
plot.data3$variable <- "normalization2"
plot.data.final <- rbind(plot.data1, rbind(plot.data2, plot.data3))

# get marker summary boxplot
jpeg(boxplot.marker.filename, width = 4, height = 5, units = "in", res = 400)
ggplot(plot.data.final, aes(x = variable, y = value)) + geom_boxplot() + ylab(halo.marker.yaxis.label)
dev.off()

# get marker summary boxplot, logscale
jpeg(boxplot.log.marker.filename)
ggplot(plot.data.final, aes(x = variable, y = value)) + geom_boxplot() + scale_y_log10() +
  ylab(paste0(halo.marker.yaxis.label, " (log10)"))
dev.off()

# histogram density plot -- nucleus area
jpeg(histogram.filename.area, width = 4, height = 2, units = "in", res = 400)
halo.outputs.table[,cell.area.variable] %>% hist() %>% plot(xlab = cell.area.variable, main = "")
dev.off()

# histogram density plot -- akt3 copies
jpeg(histogram.filename.marker, width = 4, height = 2, units = "in", res = 400)
halo.outputs.table[,gene.marker.label] %>% hist() %>% plot(xlab = gene.marker.label, main = "")
dev.off()
