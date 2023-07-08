#!/usr/bin/env R

# Author: Sean Maden
#
# Concordance at the top plots for marker sets.
#


libv <- c("ggplot2", "ggforce", "gridExtra", "dplyr")
sapply(libv, library, character.only = TRUE)

set.seed(0) # seed for random colors

#----------
# load data
#----------
# ro1 dlpfc sce, markers
fname1 <- "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda"
lscef1 <- get(load(fname1))

# mrb dlpfc sce, markers
fname2 <- "list-scef_markers-k2-k3-k4_mrb-dlpfc.rda"
lscef2 <- get(load(fname2))

#-----------------
# helper functions
#-----------------

cat_dfp <- function(subject.list, compare.list){
  # cat_dfp
  #
  # Get the plot table for a CAT plot.
  #
  # subject.list : List containing discrete marker sets for overlaps.
  # compare.list : List of individual sets to compare with subject.list.
  #
  sl <- subject.list; cl <- compare.list; clnames <- names(cl)
  dfp <- do.call(rbind, lapply(seq(length(sl)), function(ii){
    sli <- sl[[ii]]; indexi <- length(sli)
    olv <- unlist(lapply(cl, function(cli){
      length(intersect(cli, sli))
    }))
    dfoli <- data.frame(compare.name = clnames, overlaps = olv)
    dfoli$subject.index <- indexi; dfoli
  }))
  dfp <- as.data.frame(dfp)
  for(ii in c(2:3)){dfp[,ii] <- as.numeric(dfp[,ii])}
  return(dfp)
}


mri <- metadata(lscef1$k2)[["markers"]][["all"]]
# get the max markers by type as min total markers among subject labels
dff <- as.data.frame(table(mri$cellType.target))
max.markers <- min(dff[,2])
# get subject.list as sequences of top markers
typev <- dff[,1]
subject.list <- lapply(seq(max.markers), function(nmarker){
  message("Getting ",nmarker, " top markers per type...")
  mrtop <- do.call(rbind, lapply(typev, function(typei){
    mri %>% filter(cellType.target == typei) %>% 
      arrange(rank_ratio) %>% top_n(n = nmarker)
  }))
  mrtop$gene
})
names(subject.list) <- seq(max.markers)

# get compare.list
compare.list <- list()
compare.list[["k3_ro1"]] <- metadata(lscef1[["k3"]])[["markers"]][["top"]]$gene
compare.list[["k4_ro1"]] <- metadata(lscef1[["k4"]])[["markers"]][["top"]]$gene
compare.list[["k2_mrb"]] <- metadata(lscef2[["k2"]])[["markers"]][["top"]]$gene
compare.list[["k3_mrb"]] <- metadata(lscef2[["k3"]])[["markers"]][["top"]]$gene
compare.list[["k4_mrb"]] <- metadata(lscef2[["k4"]])[["markers"]][["top"]]$gene

# get dfp
dfp <- cat_dfp(subject.list, compare.list)
dfp$Comparison <- dfp$compare.name

# get plot
title.str <- "CAT: K2 RO1 DLPFC"
ggcat <- ggplot(dfp, aes(x = subject.index, y = overlaps, color = Comparison)) + 
  geom_point() + geom_line() + theme_bw() + ggtitle(title.str) +
  ylab("Overlapping markers") + xlab("Top marker count") +
  facet_zoom(ylim = c(0, 20), xlim = c(0, 50), zoom.size = 1, show.area = FALSE)
# save plot
save.dpath <- "."
fname <- "ggcat_k2-markers_ro1-dlpfc.jpg"
jpeg(file.path(save.dpath, fname), 
     width = 6, height = 2.5, units = "in", res = 400)
print(ggcat); dev.off()

