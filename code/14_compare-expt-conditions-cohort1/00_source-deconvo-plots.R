#!/usr/bin/env R

# Author: Sean Maden
#
# Deconvolution plots script. Sourced by notebook for Ybulk experiment with 4 identifier conditions in the experiment.
#

highlight.color.low = "red"
highlight.color.high = "blue"
background.color.highlights = "gray"

deconvo_plot_statistics <- function(dfp1, error.type = "type1"){
  #
  #
  #
  
  # postprocess dfp1
  
  # append fractions
  dfp1$s.fraction.type2.type1 <- dfp1$s.type2/dfp1$s.type1
  dfp1$s.fraction.type2.type1.discrete <- as.character(round(dfp1$s.fraction.type2.type1, 2))
  dfp1$log.s.fraction <- log(dfp1$s.fraction.type2.type1)
  
  # append error summaries
  err.variable.name <- paste0("error.", error.type)
  dfp1$error.value <- err.variable <- dfp1[,err.variable.name]
  dfp1$minimum.error <- err.variable==min(err.variable)
  dfp1$maximum.error <- err.variable==max(err.variable)
  deciles.error <- quantile(err.variable, seq(0, 1, 0.1))
  dfp1$minimum.decile.error <- err.variable <= deciles.error[2]
  dfp1$maximum.decile.error <- err.variable >= deciles.error[9]
  
  
  # get highlight categories
  dfp1$type1.group.label <- as.character(dfp1$s.type1)
  dfp1$type2.group.label <- as.character(dfp1$s.type2)
  dfp1$all.highlight.categories <- ifelse(dfp1$minimum.error, "min", 
                                            ifelse(dfp1$maximum.error, "max", 
                                                   ifelse(dfp1$minimum.decile.error, "min.dec", 
                                                          ifelse(dfp1$maximum.decile.error, "max.dec", "mid"))))
  dfp1$all.highlight.sizes <- ifelse(dfp1$minimum.error|dfp1$maximum.error, "min/max", "min.dec/max.dec/mid")
  
  return(dfp1)
}

deconvo_heatmaps <- function(dfp.iter, error.type = "type1", draw.min.err.line = TRUE,
                             xlab.string = "type1", ylab.string = "type2"){
  #
  #
  #
  #
  
  # make the plots
  # Plot heatmap -- error.type2
  heatmap1 <- ggplot(dfp.iter, aes(x = s.type1, y = s.type2)) + theme_bw() +
    geom_raster(aes(fill = error.value)) +
    scale_fill_gradientn(colours = rainbow(5)) + 
    geom_abline(slope = 1, intercept = 0) +
    ggtitle(dfp.iter$block.experiment) + 
    ggtitle("Error: ", error.type) +
    xlab(xlab.string) + ylab(ylab.string)
  
  # prep heatlution plots with highlights
  # recompute min/max/decile errors
  dfp.iter$minimum.error <- dfp.iter$error.value==min(dfp.iter$error.value)
  dfp.iter$maximum.error <- dfp.iter$error.value==max(dfp.iter$error.value)
  deciles.error.vector <- quantile(dfp.iter$error.value, seq(0, 1, 0.1))
  dfp.iter$minimum.decile.error <- dfp.iter$error.value <= deciles.error.vector[2]
  dfp.iter$maximum.decile.error <- dfp.iter$error.value >= deciles.error.vector[10]
  
  # min highlights
  heatmap2 <- ggplot(dfp.iter, aes(x = s.type1, y = s.type2, color = minimum.error)) + 
    geom_point(alpha = 1) + scale_color_manual(values = c("TRUE" = highlight.color.low, 
                                                          "FALSE" = background.color.highlights)) +
    ggtitle("Error: ", error.type) + xlab(xlab.string) + ylab(ylab.string)
  heatmap3 <- ggplot(dfp.iter, aes(x = s.type1, y = s.type2, color = minimum.decile.error)) + 
    geom_point() + scale_color_manual(values = c("TRUE" = highlight.color.low, 
                                                 "FALSE" = background.color.highlights)) +
    ggtitle("Error: ", error.type) + xlab(xlab.string) + ylab(ylab.string)
  # max highlights
  heatmap4 <- ggplot(dfp.iter, aes(x = s.type1, y = s.type2, color = maximum.error)) + 
    geom_point() + scale_color_manual(values = c("TRUE" = highlight.color.high, 
                                                 "FALSE" = background.color.highlights)) +
    ggtitle("Error: ", error.type) + xlab(xlab.string) + ylab(ylab.string)
  heatmap5 <- ggplot(dfp.iter, aes(x = s.type1, y = s.type2, color = maximum.decile.error)) + 
    geom_point() + scale_color_manual(values = c("TRUE" = highlight.color.high, 
                                                 "FALSE" = background.color.highlights)) +
    ggtitle("Error: ", error.type) + xlab(xlab.string) + ylab(ylab.string)
  heatmap6 <- ggplot(dfp.iter, 
                     aes(x = s.type1, y = s.type2, color = all.highlight.categories, 
                               size = all.highlight.sizes)) + 
    geom_abline(slope = 1, intercept = 0) +
    geom_point() + 
    scale_color_manual(values = 
                         c("mid" = background.color.highlights, 
                           "min" = highlight.color.low,
                           "max" = highlight.color.high,
                           "min.dec" = highlight.color.low,
                           "max.dec" = highlight.color.high)) +
    scale_size_manual(values = c("min/max" = 2.5, 
                                 "min.dec/max.dec/mid" = 1)) +
    ggtitle("Error: ", error.type) + xlab(xlab.string) + ylab(ylab.string)
  
  # parse min.err line
  if(draw.min.err.line){
    # get plot lines
    dfp.iter.min.err <- dfp.iter[dfp.iter$minimum.error,]
    slope.min.err <- dfp.iter.min.err$s.type2/dfp.iter.min.err$s.type1
    # get revised plots
    heatmap1 <- heatmap1 + geom_abline(slope = slope.min.err, intercept = 0, linetype = "dashed", lwd = 2)
    heatmap2 <- heatmap2 + geom_abline(slope = slope.min.err, intercept = 0, linetype = "dashed", lwd = 2)
    heatmap3 <- heatmap3 + geom_abline(slope = slope.min.err, intercept = 0, linetype = "dashed", lwd = 2)
    heatmap4 <- heatmap4 + geom_abline(slope = slope.min.err, intercept = 0, linetype = "dashed", lwd = 2)
    heatmap5 <- heatmap5 + geom_abline(slope = slope.min.err, intercept = 0, linetype = "dashed", lwd = 2)
    heatmap6 <- heatmap6 + geom_abline(slope = slope.min.err, intercept = 0, linetype = "dashed", lwd = 2)
  }

  # get return list
  list(heatmap1 = heatmap1, heatmap2 = heatmap2, heatmap3 = heatmap3, 
       heatmap4 = heatmap4, heatmap5 = heatmap5, heatmap6 = heatmap6)
}

deconvo_scatterplots <- function(dfp.iter){
  # SCATTERLUTION PLOTS
  sp1 <- ggplot(dfp.iter, aes(x = s.type2, y = bias.type2, group = s.type1)) + 
    geom_point(alpha = 1) + geom_line() + theme_bw() + 
    facet_wrap(~sample.id) + geom_hline(yintercept = 0)
  # Plot points with s.fraction value as color
  # continuous grouping
  sp2 <- ggplot(dfp.iter, aes(x = s.type2, y = bias.type2, 
                          color = s.fraction.type2.type1, 
                          group = s.fraction.type2.type1)) + 
    geom_point(size = 2) + theme_bw() + 
    facet_wrap(~sample.id) + geom_hline(yintercept = 0)
  list(sp1 = sp1, sp2 = sp2)
}

deconvo_lineplots <- function(dfp.iter){
  # DECONVOLINE PLOTS
  lp1 <- ggplot(dfp.iter, aes(x = s.type2, y = bias.type2.true.pred, group = s.type1)) + 
    geom_line() + theme_bw() + 
    facet_wrap(~sample.id) + geom_hline(yintercept = 0)
  lp2 <- ggplot(dfp.iter, aes(x = s.type2, y = bias.type2.true.pred, group = s.type1, color = s.type1)) + 
    geom_line() + theme_bw() + 
    facet_wrap(~sample.id) + geom_hline(yintercept = 0)
  # discrete grouping
  lp3 <- ggplot(dfp.iter, aes(x = s.type2, y = bias.type2.true.pred, 
                              group = s.type1.group.label, color = s.type1.group.label)) + 
    geom_line() + theme_bw() + 
    facet_wrap(~sample.id) + geom_hline(yintercept = 0)
  lp4 <- ggplot(dfp.iter, aes(x = s.type2, y = bias.type2.true.pred, 
                              group = s.type1.group.label, color = s.type1.group.label)) + 
    geom_line() + theme_bw() + 
    facet_wrap(~sample.id) + geom_hline(yintercept = 0)
  list(lp1 = lp1, lp2 = lp2, lp3 = lp3, lp4 = lp4)
}

deconvocanos <- function(dfp.iter){
  # DECONVOCANO_PLOTS
  # Scatterplot of fraction versus bias.
  cano1 <- ggplot(dfp.iter, aes(x = log.s.fraction, y = bias.type2.true.pred)) + 
    geom_point(size = 2) + theme_bw() + facet_wrap(~sample.id) + geom_hline(yintercept = 0)
  # Scatterplot of fraction versus error
  cano2 <- ggplot(dfp.iter, aes(x = log.s.fraction, y = error.type2)) + 
    geom_point(size = 2) + theme_bw() + facet_wrap(~sample.id) + geom_hline(yintercept = 0)
  list(cano1 = cano1, cano2 = cano2, metadata = list(plottype="volcano_plot"))
}

#
# viz wrapper function -- NEEDS WORK
# 
deconvo_plots_list <- function(dfp, facet.variable = NULL, 
                               update.stat.summaries = TRUE,
                               draw.min.err.line = TRUE){
  # get cell types vector
  colnames.dfp <- colnames(dfp)
  colnames.types <- colnames.dfp[1:(which(grepl("sample.label", colnames.dfp))-1)]
  types.vector <- gsub("\\..*", "", colnames.types)
  types.combos <- combn(types.vector, 2)
  
  # get plots of pairwise cell type comparisons
  plots.list <- lapply(seq(ncol(types.combos)), function(index){
    message(index)
    combo.iter <- types.combos[,index]
    type1.iter <- combo.iter[1]
    type2.iter <- combo.iter[2]
    # get plot data
    s.type1.variable <- dfp[,grepl(paste0("^s.", type1.iter, ".*"), colnames.dfp)]
    s.type2.variable <- dfp[,grepl(paste0("^s.", type2.iter, ".*"), colnames.dfp)]
    bias.type1.variable <- dfp[,grepl(paste0("^bias.", type1.iter, "\\..*"), colnames.dfp)]
    bias.type2.variable <- dfp[,grepl(paste0("^bias.", type2.iter, "\\..*"), colnames.dfp)]
    error.type1.variable <- dfp[,grepl(paste0("^error.", type1.iter, "\\..*"), colnames.dfp)]
    error.type2.variable <- dfp[,grepl(paste0("^error.", type2.iter, "\\..*"), colnames.dfp)]
    dfp.iter <- data.frame(s.type1 = s.type1.variable,
                           s.type2 = s.type2.variable,
                           bias.type1.true.pred = bias.type1.variable,
                           bias.type2.true.pred = bias.type2.variable,
                           error.type1 = error.type1.variable,
                           error.type2 = error.type2.variable,
                           sampe.label = dfp$sample.label,
                           sample.id = dfp$sample.id)
    # calculate plot statistics
    dfp.iter <- deconvo_plot_statistics(dfp.iter, error.type = "type1")
    # gets the plots
    lr <- list(heatmaps = 
                 deconvo_heatmaps(dfp.iter, type1.iter, draw.min.err.line = T,
                                  xlab.string = type1.iter, ylab.string = type2.iter),
               scatterplots = deconvo_scatterplots(dfp.iter),
               lineplots = deconvo_lineplots(dfp.iter),
               deconvocanos = deconvocanos(dfp.iter),
               dfp = dfp.iter)
    # parse facet option on ALL OF THE PLOTS
    if(is(facet.variable, "NULL")){
    } else{
      lr.labels <- names(lr)
      lr <- lapply(lr, function(lr.plot.type){
        plot.labels <- names(lr.plot.type)
        lr.plot.type <- lapply(lr.plot.type, function(plot){
          plot <- eval(parse(text = paste0("plot + facet_wrap(~",facet.variable,")")))
        })
        names(lr.plot.type) <- plot.labels
        return(lr.plot.type)
      })
      names(lr) <- lr.labels
    }
    lr$metadata <- list(type.label.combo = combo.iter)
    return(lr)
  })
  names(plots.list) <- unlist(apply(types.combos,2,paste0,collapse = ";"))
  # return
  lr <- plots.list
  return(lr)
}

unpack_plot_type <- function(deconvo_plots_list_list, plot.name, plot.type){
  # this is for lists of deconvo_plots_list objects
  #
  # notes:
  # * replaces: lapply(results.list.dfp4, function(results.list.iter){results.list.iter$heatmaps$heatmap1})
  #
  lapply(results.list.dfp4, function(results.list.iter){results.list.iter[[plot.type]][[plot.name]]})
}