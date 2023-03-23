#!/usr/bin/env R

# Author: Sean Maden
#
# Run donor bias subsampling experiments.
#
# See also: script 05
#

source("deconvo_method-paper/code/10_donor-bias-simulations-continued/00_parameters.R")
sapply(libv, library, character.only = T)

# analyze results table
# get results table
data.directory.files <- list.files(base.path.workflow.data)
filter.files <- grepl(results.filt, list.data.directory)
results.table.file.name <- data.directory.files[filter.files]
results.table.path <- file.path(base.path.workflow.data, results.table.file.name)
results.table <- read.csv(results.table.path)
# remove duplicate entries
results.table.filtered <- results.table
results.table.filtered$label <- paste0(results.table.filtered[,iterations.index.column.name], 
                                       ";", 
                                       results.table.filtered[,algorithm.column.name])
results.table.filtered <- results.table.filtered[!duplicated(results.table.filtered$label),]
methodv <- unique(rtf[,method.varname])

# make new bias plots
# get list of jitter plots
ggplot.list <- lapply(variable.vector, function(variable){
  plot.data <- results.table.filtered
  plot.data$value <- plot.data[,variable]
  plot.data$x.variable <- plot.data[,algorithm.column.name]
  ggplot(plot.data, aes(x = x.variable, y = value)) +
    geom_jitter(alpha = 0.5) + ggtitle(variable) +
    geom_boxplot(alpha = 0, color = "cyan", size = 1) +
    xlab(jitter.plot.xaxis.label)
})
# get plots from list
plot1 <- ggplot.list[[3]] + theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  geom_hline(yintercept = 0, alpha = 1, color = "red")
plot2 <- ggplot.list[[4]] + theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, alpha = 1, color = "red")
# save
jpeg("ggjitterbox-comp-bias_inter-sample-bias.jpg", 
     width = 6, height = 5, units = "in", res = 400)
grid.arrange(plot1, plot2, nrow = 2,
             layout_matrix = matrix(c(1,1,1,2,2,2,2), ncol = 1),
             bottom = "Deconvolution method")
dev.off()

# plot pairs
# varname <- "prop.pred.type1"
# iterate on plot variables
ggplot.list <- lapply(variable.vector, function(variable.name){
  plot.data.frame <- do.call(cbind, lapply(algorithm.vector, function(algorithm){
    filter.results <- rtf[,algorithm.column.name]==algorithm
    results.table.filtered.iter <- results.table.filtered[filter.results,]
    results.table.filtered.iter <- results.table.filtered.iter[
      order(results.table.filtered.iter$iterations_index),]
    results.table.filtered.iter[,variable.name]
  }))
  plot.data.frame <- plot.data.frame %>% as.data.frame()
  colnames(plot.data.frame) <- algorithm.vector
  ggpairs(plot.data.frame) + ggtitle(variable.name)
})
names(ggplot.list) <- variable.vector
# show new plots
ggplot.list$prop.pred.type1
ggplot.list$prop.pred.type2
ggplot.list$bias.type1
ggplot.list$bias.type2
ggplot.list$rmse.types

# get five violin plots
ggplot.list <- lapply(varv, function(varname){
  dfp <- rtf; dfp$value <- dfp[,varname]
  ggplot(dfp, aes(x = deconvolution_method, y = value)) +
    geom_violin(draw_quantiles = 0.5) + ggtitle(varname)
})
jpeg(violin.plots.five.filename, width = 5, height = 4, units = "in", res = 400)
grid.arrange(ggplot.list[[1]], ggplot.list[[2]], ggplot.list[[3]], ggplot.list[[4]], ggplot.list[[5]], layout_matrix = matrix(c(1,2,3,4,5), nrow = 2))
dev.off()

# get two violin plots
plot1 <- ggplot.list[[3]] + theme_bw() + theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + geom_hline(yintercept = 0, alpha = 0.2)
plot2 <- ggplot.list[[4]] + theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + geom_hline(yintercept = 0, alpha = 0.2)
jpeg(violin.plots.two.filename, width = 6, height = 5, units = "in", res = 400)
grid.arrange(plot1, plot2, nrow = 2, layout_matrix = matrix(c(1,1,1,2,2,2,2), ncol = 1), bottom = "Deconvolution method")
dev.off()

# scatter plot bias
# no facet
ggpt <- ggplot(results.table.filtered, aes(x = bias.type1, y = bias.type2, group = deconvolution_method, color = deconvolution_method)) +
  geom_point(alpha = 0.2) + 
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() + 
  geom_hline(yintercept = 0, alpha = 0.5, color = "gray") +
  geom_vline(xintercept = 0, alpha = 0.5, color = "gray")
# facet
ggpt <- ggplot(rtf, aes(x = bias.type1, y = bias.type2)) +
  geom_point(alpha = 0.2) + 
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() + 
  geom_hline(yintercept = 0, alpha = 0.5, color = "gray") +
  geom_vline(xintercept = 0, alpha = 0.5, color = "gray")
ggpt <- ggpt + facet_wrap(~deconvolution_method)
jpeg(scatterplot.filename, width = 6, height = 6, units = "in", res = 400)
ggpt
dev.off()
