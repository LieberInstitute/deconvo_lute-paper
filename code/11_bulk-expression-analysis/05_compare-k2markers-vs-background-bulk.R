#!/usr/bin/env R

# Author: Sean Maden
#
#

# load data
source("deconvo_method-paper/code/11_bulk-expression-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
rse.background <- get(load(rse.gene.filter.filepath))
rse.markers <- get(load(rse.k2markers.filepath))
cd <- colData(rse.background)

# prepare comparison params
# set up data summaries
# params

# get expression scales to analyze
list.comparison.plots <- lapply(assays, function(assay.name){
  expression.background <- assays(rse.background)[[assay.name]]
  expression.markers <- assays(rse.markers)[[assay.name]]
  get_comparison_data(expression.background = expression.background, 
                      expression.markers = expression.markers,
                      plot.filename = plot.filename, 
                      variable.vector = variable.vector,
                      type.vector = type.vector,
                      cd = cd, save.path = save.path)
})

# save new plots
# get logscale boxplots list
plot.list <- list.comparison.plots[[1]]$lgg$box.log
statistic.type.list <- names(plot.list)
statistic.index.vector <- seq(length(statistic.type.list))
# save composite plots
for(index in statistic.index.vector){
  statistic.name <- statistic.type.list[index]
  plot.path <- bulk.marker.compare.plots.paths[index]
  plot.list.iter <- plot.list[[index]]
  
  # get formatted plots by group type
  plot1 <- plot.list.iter[[1]] + 
    scale_colour_manual(values = colors.compare.markers) +
    theme(axis.title.y = element_blank(), 
          axis.title.x = element_blank())
  plot2 <- plot.list.iter[[2]] + 
    scale_colour_manual(values = colors.compare.markers) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank())
  plot3 <- plot.list.iter[[3]] + 
    scale_colour_manual(values = colors.compare.markers) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank())
  plot4 <- plot.list.iter[[4]] + 
    scale_colour_manual(values = colors.compare.markers) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank())
  
  # save
  jpeg(plot.path, width = 10, height = 5, units = "in", res = 400)
  grid.arrange(plot1, plot2, plot3, plot4, nrow = 2, left = statistic.name)
  dev.off()
}





# qc at markers versus background -- jitter/boxplots
# mean expression
# expression variance
# correlation bulk and halo marker expression
# marker degs among bulk samples