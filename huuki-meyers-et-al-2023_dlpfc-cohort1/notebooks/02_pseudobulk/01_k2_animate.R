#!/usr/bin/env R

# Author: Sean Maden
#
# Makes gganimate animation from pseudobulk results with dot-and-line segement for
# two cell types.
#
#

libv <- c("ggplot2", "gganimate", "ggtext")
sapply(libv, library, character.only = TRUE)
knitr::opts_chunk$set(echo = TRUE)
list.files()
load("./env/02_pseudobulk/01_k2.RData")

# make plot data
dfi <- dfp.tall
dfi1 <- data.frame(
  true = dfi$neuron.true,
  pred = dfi$neuron.pred,
  sample.id = dfi$sample.id,
  type = dfi$type,
  celltype = "neuron"
)
dfi2 <- data.frame(
  true = dfi$glial.true,
  pred = dfi$glial.pred,
  sample.id = dfi$sample.id,
  type = dfi$type,
  celltype = "glial"
)
dfp <- rbind(dfi1, dfi2) |> as.data.frame()

# make plots
cellTypeColor1 <- "blue"
cellTypeColor2 <- "red"
idLineColor <- "orange"

# get plot objects
new.plot <- ggplot(dfp, aes(x = true, y = pred, colour = celltype)) + 
  geom_point(size = 4, alpha = 0.3) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_hline(yintercept = 0.5) + 
  geom_vline(xintercept = 0.5) + theme_bw() +
  xlab("Known") + ylab("Predicted") + xlim(0, 1) + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c(cellTypeColor2, cellTypeColor1))
new.plot <- new.plot + geom_line(aes(group=sample.id), 
                                 color = idLineColor, alpha = 0.5)
new.plot

# make animation
ani <- new.plot + 
  labs(title = "DLPFC pseudobulk ({closest_state})") + 
  transition_states(type)

animate(
  ani,
  width = 3.5,
  height = 2.8,
  units = "in",
  res = 400,
  renderer = gifski_renderer(file = "fig2_dlpfc-k2.gif")
)
