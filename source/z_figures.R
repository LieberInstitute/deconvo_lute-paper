#!/usr/bin/env R

#
# Functions to make standard figures for z tables (e.g. zsource, zfinal, z_dk, 
# etc.).
#

get_lgg_markers <- function(df.markers, save.new.plots = T,
                            save.dpath = "deconvo_method-paper/outputs/02_test-source/",
                            png.fname.cplot = "mr-k_composite-vp-jt-bp.png",
                            plot.single.fnamev = c("mr-k_violin-plot.png", 
                                                   "mr-k_jitter-plot.png",
                                                   "mr-k_box-plot.png"),
                            png.width.splot = 4,
                            png.height.splot = 3,
                            png.units = "in",
                            png.res = 400,
                            png.width.cplot = 4,
                            png.height.cplot = 6){
  # Get a list of ggplots of gene marker mean ratio summaries. Also makes new
  # figures.
  #
  # df.markers: "top.marker.data" object from lz outputs. This should have been
  #   output from DeconvoBuddies::get_mean_ratio2().
  # save.new.plots: Whether to save new plots of ggplot objects.
  #
  #
  dfpi <- as.data.frame(df.markers[,c(1,2,4,6,7)]) # index colnames of interest
  # mean ratio by cell type
  # plots for individual figures
  message("making objects for individual plots...")
  ggvp1 <- ggplot(dfpi, aes(x = celltype.target, y = ratio, fill = celltype.target)) + 
    theme_bw() + geom_violin(draw_quantiles = 0.5) + ylab("MR_k") + xlab("cell_type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  ggjt1 <- ggplot(dfpi, aes(x = celltype.target, y = ratio, color = celltype.target)) + 
    theme_bw() + ylab("MR_k") + xlab("cell_type") + geom_jitter(size = 1, alpha = 0.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  ggbp1 <- ggplot(dfpi, aes(x = celltype.target, y = ratio, fill = celltype.target)) + 
    theme_bw() + geom_boxplot() + ylab("MR_k") + xlab("cell_type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  # plots for composite figures
  message("making objects for composite plot...")
  ggvp2 <- ggplot(dfpi, aes(x = celltype.target, y = ratio, fill = celltype.target)) + 
    theme_bw() + geom_violin(draw_quantiles = 0.5) + ylab("MR_k") + xlab("cell_type") +
    theme(axis.text.x = element_blank(), legend.position = "none",
          axis.title.x = element_blank(), axis.title.y = element_blank()) +
    ggtitle("Violin plots")
  ggjt2 <- ggplot(dfpi, aes(x = celltype.target, y = ratio, color = celltype.target)) + 
    theme_bw() + ylab("MR_k") + xlab("cell_type") + geom_jitter(size = 1, alpha = 0.5) +
    theme(axis.text.x = element_blank(), legend.position = "none",
          axis.title.x = element_blank(), axis.title.y = element_blank()) +
    ggtitle("Jitter points")
  ggbp2 <- ggplot(dfpi, aes(x = celltype.target, y = ratio, fill = celltype.target)) + 
    theme_bw() + geom_boxplot() + ylab("MR_k") + xlab("cell_type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          legend.position = "none") +
    ggtitle("Boxplots")
  # save plots
  if(save.new.plots){
    message("Saving new plot files...")
    # save individual plots
    plot.single.fnamev <- c("mr-k_violin-plot.png", "mr-k_jitter-plot.png",
                            "mr-k_box-plot.png")
    plot.fpathv <- file.path(save.dpath, plot.single.fnamev)
    for(ii in seq(length(lr[[1]]))){
      message("saving new plot: ", plot.fpathv[ii])
      png(plot.fpathv[ii], width = png.width.splot, 
          height = png.height.splot, res = png.res, units = png.units)
      print(lr[[1]][[ii]]); dev.off()
    }
    # save composite plot
    cplot.fpath <- file.path(save.dpath, png.fname.cplot)
    message("saving new plot: ", cplot.fpath)
    lm <- matrix(c(1,1,2,2,3,3,3),ncol=1) # adds extra space for xaxis labels
    png(cplot.fpath, width = png.width.cplot, 
        height = png.height.cplot, units = png.units, res = png.res)
    grid.arrange(ggvp2, ggjt2, ggbp2, ncol = 1, layout_matrix = lm, 
                 left = paste0(paste0(rep(" ", 20),collapse=""), "MR_k"), 
                 bottom = "cell_type")
    dev.off()
  }
  # get return list
  lr <- list(single.plots = list(violin = ggvp1, jitter = ggjt1, box = ggbp1),
             composite.plots = list(violin = ggvp2, jitter = ggjt2, box = ggbp2))
  return(lr)
}