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

# make summary  plots
# violin plot
ggplot(plot.data.final, aes(x = variable, y = value)) + geom_violin(draw_quantiles = 0.5)
# boxplot
ggplot(plot.data.final, aes(x = variable, y = value)) + geom_boxplot()
# boxplot with logscale
ggplot(plot.data.final, aes(x = variable, y = value)) + geom_boxplot() + scale_y_log10()

# check normality
# plot histograms
halo.outputs.table[,cell.area.variable] %>% hist() %>% plot()
halo.outputs.table[,gene.marker.label] %>% hist() %>% plot()
vector <- max(halo.outputs.table[,gene.marker.label] + 1) - halo.outputs.table[,gene.marker.label]
transformed.halo.outputs <- 1/vector
halo.outputs.table[,gene.marker.label] %>% log() %>% hist() %>% plot()

