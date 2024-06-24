#!/usr/bin/env R

# Author: Sean Maden
#
# Makes gganimate animation from real bulk sample NNLS deconvolution results 
# with dot-and-line segement for two cell types.
#
#

knitr::opts_chunk$set(echo = TRUE)
libv <- c("dplyr", "ggplot2", "GGally", "gridExtra", "gganimate", "ggtext")
sapply(libv, library, character.only = TRUE)
load("./env/07_adjustment/03_run_adjustment_realbulk_all_script.RData")

# get plot data
# parse dfp.tall
dfp.tall <- dfp.tall[dfp.tall$cell.type=="neuron",]
# append glial data
# pred value
dfp.tall2 <- dfp.tall
dfp.tall2$cell.type <- "glial"
dfp.tall2$value <- 1-dfp.tall2$value
dfp.tall2$true <- 1-dfp.tall2$true
# error
dfp.tall2$error.negative <- dfp.tall2$value-dfp.tall2$true
dfp.tall2$error <- abs(dfp.tall2$error.negative)
# recombine
# bind
dfp.tall$error.negative <- dfp.tall$true-dfp.tall$value
dfp.tall <- rbind(dfp.tall, dfp.tall2) |> as.data.frame()
# rename type variables
dfp.tall$type <- ifelse(dfp.tall$scale, "withscale", "noscale")

# make plot objects
cellTypeColor1 <- "blue"
cellTypeColor2 <- "red"
idLineColor <- "orange"
new.plot <- ggplot(dfp.tall[dfp.tall$algorithm=="nnls",], 
                                  aes(x = true, y = value, colour = cell.type)) + 
  geom_point(size = 4, alpha = 0.3) + geom_abline(slope = 1, intercept = 0) + 
  geom_hline(yintercept = 0.5) + geom_vline(xintercept = 0.5) + theme_bw() +
  xlab("Known") + ylab("Predicted") + xlim(0, 1) + ylim(0, 1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c(cellTypeColor2, cellTypeColor1))
new.plot <- new.plot + geom_line(aes(group=sample.id), color = idLineColor,
                                 alpha = 0.5)
new.plot

# make and save new animation
ani <- new.plot + 
  labs(title = "DLPFC bulk ({closest_state})") + 
  transition_states(type)
animate(
  ani,
  width = 3.5,
  height = 3,
  units = "in",
  res = 400,
  renderer = gifski_renderer(
    file = "fig3_dlpfc-realbulk-k2.gif"
  )
)