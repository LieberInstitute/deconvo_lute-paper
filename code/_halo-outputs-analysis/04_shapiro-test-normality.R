#!/usr/bin/env R

# Author: Sean Maden
#

source("deconvo_method-paper/code/_halo-outputs-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
halo.outputs.table <- get(load(halo.output.path))
halo.outputs.table <- halo.outputs.table %>% as.data.frame()

# shapiro-wilk tests of normality
test.input <- gene.marker.vector
shapiro.test.transform1 <- sample(test.input, shapiro.downsample.amount) %>%shapiro.test()

test.input <- 1/(halo.outputs.table[,gene.marker.label]+1e-10)
shapiro.test.transform2 <- sample(test.input, shapiro.downsample.amount) %>% shapiro.test()

test.input <- sqrt(gene.marker.vector)
shapiro.test.transform3 <- sample(test.input, shapiro.downsample.amount) %>% shapiro.test()
