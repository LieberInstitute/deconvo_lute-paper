#!/usr/bin/env R

# Author: Sean Maden
#
# Deconvolution plots script. Sourced by notebook for Ybulk experiment with 4 identifier conditions in the experiment.
#

deconvo_plot_statistics <- function(dfp1){
  # postprocess dfp1
  dfp1$s.fraction.neuron.glial <- dfp1$neuron/dfp1$glial
  dfp1$s.fraction.neuron.glial.discrete <- as.character(round(dfp1$s.fraction.neuron.glial, 2))
  dfp1$log.s.fraction <- log(dfp1$s.fraction.neuron.glial)
  dfp1$minimum.error <- dfp1$error.neuron==min(dfp1$error.neuron)
  dfp1$maximum.error <- dfp1$error.neuron==max(dfp1$error.neuron)
  deciles.error.neuron <- quantile(dfp1$error.neuron, seq(0, 1, 0.1))
  dfp1$minimum.decile.error <- dfp1$error.neuron <= deciles.error.neuron[2]
  dfp1$maximum.decile.error <- dfp1$error.neuron >= deciles.error.neuron[9]
  dfp1$glial.group.label <- as.character(dfp1$glial)
  dfp1$neuron.group.label <- as.character(dfp1$neuron)
  # dfp1$all.highlight.categories <- ifelse(dfp1$minimum.error, "min",ifelse(,,ifelse(,,ifelse())))
  dfp1$all.highlight.categories <- ifelse(dfp1$minimum.error, "min", 
                                            ifelse(dfp1$maximum.error, "max", 
                                                   ifelse(dfp1$minimum.decile.error, "min.dec", 
                                                          ifelse(dfp1$maximum.decile.error, "max.dec", "mid"))))
  dfp1$all.highlight.sizes <- ifelse(dfp1$minimum.error|dfp1$maximum.error, "min/max", "min.dec/max.dec/mid")
  return(dfp1)
}

deconvo_heatmaps <- function(dfp1, facet.variable = NULL){
  #---------------
  # make the plots
  #---------------
  # Plot heatmap -- error.neuron
  heatmap1 <- ggplot(dfp1, aes(x = glial, y = neuron)) + theme_bw() +
    geom_raster(aes(fill = error.neuron)) +
    scale_fill_gradientn(colours = rainbow(5)) + 
    geom_vline(xintercept = 3) + geom_hline(yintercept = 10) +
    ggtitle(dfp1$block.experiment)
  
  # prep heatlution plots with highlights
  # recompute min/max/decile errors
  dfp1$minimum.error <- dfp1$error.neuron==min(dfp1$error.neuron)
  dfp1$maximum.error <- dfp1$error.neuron==max(dfp1$error.neuron)
  deciles.error.neuron <- quantile(dfp1$error.neuron, seq(0, 1, 0.1))
  dfp1$minimum.decile.error <- dfp1$error.neuron <= deciles.error.neuron[2]
  dfp1$maximum.decile.error <- dfp1$error.neuron >= deciles.error.neuron[10]
  # min highlights
  heatmap2 <- ggplot(dfp1, aes(x = glial, y = neuron, color = minimum.error)) + 
    geom_point(alpha = 1) + scale_color_manual(values = c("TRUE" = highlight.color.low, 
                                                          "FALSE" = background.color.highlights))
  heatmap3 <- ggplot(dfp1, aes(x = glial, y = neuron, color = minimum.decile.error)) + 
    geom_point() + scale_color_manual(values = c("TRUE" = highlight.color.low, 
                                                 "FALSE" = background.color.highlights))
  # max highlights
  heatmap4 <- ggplot(dfp1, aes(x = glial, y = neuron, color = maximum.error)) + 
    geom_point() + scale_color_manual(values = c("TRUE" = highlight.color.high, 
                                                 "FALSE" = background.color.highlights))
  heatmap5 <- ggplot(dfp1, aes(x = glial, y = neuron, color = maximum.decile.error)) + 
    geom_point() + scale_color_manual(values = c("TRUE" = highlight.color.high, 
                                                 "FALSE" = background.color.highlights))
  heatmap6 <- ggplot(dfp1, 
                     aes(x = glial, y = neuron, color = all.highlight.categories, 
                               size = all.highlight.sizes)) + 
    geom_point() + 
    scale_color_manual(values = 
                         c("mid" = background.color.highlights, 
                           "min" = highlight.color.low,
                           "max" = highlight.color.high,
                           "min.dec" = highlight.color.low,
                           "max.dec" = highlight.color.high)) +
    scale_size_manual(values = c("min/max" = 2.5, 
                                 "min.dec/max.dec/mid" = 1))
  
  #----------------
  # get return list
  #----------------
  lr <- list(heatmap1 = heatmap1, heatmap2 = heatmap2, heatmap3 = heatmap3, 
             heatmap4 = heatmap4, heatmap5 = heatmap5, heatmap6 = heatmap6)
  #---------------------------------------
  # parse facet option on ALL OF THE PLOTS
  #---------------------------------------
  if(!is(lr, "NULL")){
    lr <- lapply(lr, function(plot){
      plot <- eval(parse(text = paste0("plot + facet_wrap(~",facet.variable,")")))
    })
  }
  return(lr) # return
}

deconvo_scatterplots <- function(dfp1){
  # SCATTERLUTION PLOTS
  sp1 <- ggplot(dfp1, aes(x = neuron, y = bias.neuron.true.pred, group = glial)) + 
    geom_point(alpha = 1) + geom_line() + theme_bw() + 
    facet_wrap(~sample.id) + geom_hline(yintercept = 0)
  # Plot points with s.fraction value as color
  # continuous grouping
  sp2 <- ggplot(dfp1, aes(x = neuron, y = bias.neuron.true.pred, 
                          color = s.fraction.neuron.glial, group = s.fraction.neuron.glial)) + 
    geom_point(size = 2) + theme_bw() + facet_wrap(~sample.id) + geom_hline(yintercept = 0)
  list(sp1 = sp1, sp2 = sp2)
}

deconvo_lineplots <- function(dfp1){
  # DECONVOLINE PLOTS
  lp1 <- ggplot(dfp1, aes(x = neuron, y = bias.neuron.true.pred, group = glial)) + 
    geom_line() + theme_bw() + 
    facet_wrap(~sample.id) + geom_hline(yintercept = 0)
  lp2 <- ggplot(dfp1, aes(x = neuron, y = bias.neuron.true.pred, group = glial, color = glial)) + 
    geom_line() + theme_bw() + 
    facet_wrap(~sample.id) + geom_hline(yintercept = 0)
  # discrete grouping
  lp3 <- ggplot(dfp1, aes(x = neuron, y = bias.neuron.true.pred, group = glial.group.label, color = glial.group.label)) + 
    geom_line() + theme_bw() + 
    facet_wrap(~sample.id) + geom_hline(yintercept = 0)
  lp4 <- ggplot(dfp1, aes(x = neuron, y = bias.neuron.true.pred, group = glial.group.label, color = glial.group.label)) + 
    geom_line() + theme_bw() + 
    facet_wrap(~sample.id) + geom_hline(yintercept = 0)
  list(lp1 = lp1, lp2 = lp2, lp3 = lp3, lp4 = lp4)
}

deconvocanos <- function(dfp1){
  # DECONVOCANO_PLOTS
  # Scatterplot of fraction versus bias.
  cano1 <- ggplot(dfp1, aes(x = log.s.fraction, y = bias.neuron.true.pred)) + 
    geom_point(size = 2) + theme_bw() + facet_wrap(~sample.id) + geom_hline(yintercept = 0)
  # Scatterplot of fraction versus error
  cano2 <- ggplot(dfp1, aes(x = log.s.fraction, y = error.neuron)) + 
    geom_point(size = 2) + theme_bw() + facet_wrap(~sample.id) + geom_hline(yintercept = 0)
  list(cano1 = cano1, cano2 = cano2, metadata = "GO 'CANO's!")
}

#
# viz wrapper function -- NEEDS WORK
# 
deconvo_plots_list <- function(dfp, facet.variable, script.path = "deconvo_plots.R", update.stat.summaries = TRUE){
  source(script.path)
  
  # gets the statistics summaries for your filtered data
  if(update.stat.summaries){dfp <- deconvo_plot_statistics(dfp)}
  
  # gets the plots
  list(heatmaps = deconvo_heatmaps(dfp, facet.variable),
       scatterplots = deconvo_scatterplots(dfp),
       lineplots = deconvo_lineplots(dfp),
       deconvocanos = deconvocanos(dfp))
}