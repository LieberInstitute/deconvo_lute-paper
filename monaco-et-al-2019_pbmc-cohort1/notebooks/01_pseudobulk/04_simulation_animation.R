#!/usr/bin/env R

# Author: Sean Maden
#
# Makes gganimate animation from PBMC pseudobulk results in two cell types.
#
#

libv <- c("lute", "ggplot2", "gridExtra", "reshape2", "ggforce", "gganimate", "ggtext")
sapply(libv, library, character.only = TRUE)
load("./env/01_pseudobulk/01_read_script.RData")
load("./env/01_pseudobulk/04_simulation.RData")
knitr::opts_chunk$set(echo = FALSE)

# get plot data
dfi <- dfPseudobulkA
dfi1 <- data.frame(
  true = dfi$true.plasmablasts,
  pred = dfi$Plasmablasts,
  sample.id = dfi$sample.id,
  type = dfi$condition,
  celltype = "plasmablast"
)
dfi2 <- data.frame(
  true = dfi$true.not.plasmablasts,
  pred = dfi$'Non-plasmablast',
  sample.id = dfi$sample.id,
  type = dfi$condition,
  celltype = "non-plasmablast"
)
dfp <- rbind(dfi1, dfi2) |> as.data.frame()

# get plot objects
cellTypeColor1 <- "blue"
cellTypeColor2 <- "red"
idLineColor <- "orange"

# get plot objects

new.plot <- ggplot(dfp, aes(x = true, y = pred, colour = celltype)) + 
  geom_point(size = 4, alpha = 0.3) + geom_abline(slope = 1, intercept = 0) + 
  geom_hline(yintercept = 0.5) + geom_vline(xintercept = 0.5) + theme_bw() +
  xlab("Known") + ylab("Predicted") + xlim(0, 1) + ylim(0, 1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c(cellTypeColor2, cellTypeColor1)) +
  facet_zoom(ylim = c(0, 0.15), xlim = c(0, 0.15), split = FALSE,
             show.area = FALSE) +
  geom_line(aes(group = sample.id), col=idLineColor, alpha = 0.5)
new.plot

# make and save new plot animations
ani <- new.plot + 
  labs(title = "PBMC pseudobulk ({closest_state})") + 
  transition_states(type) +
  theme(plot.title = element_markdown())
  
animate(
  ani, width = 5.5, height = 2.5, units = "in", res = 400,
  renderer = gifski_renderer(
    file = "fig2_pbmc-k2.gif"
  )
)