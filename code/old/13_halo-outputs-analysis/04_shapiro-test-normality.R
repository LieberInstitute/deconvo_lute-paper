#!/usr/bin/env R

# Author: Sean Maden
#
# Evaluate normality of log10 transformed nucleus area data.
#

source("deconvo_method-paper/code/13_halo-outputs-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
halo.outputs.table <- get(load(output.updated.path))
halo.outputs.table <- halo.outputs.table %>% as.data.frame()

# test normality
# untransformed area
test.input <- halo.outputs.table[,cell.area.variable]
shapiro.test.untransformed <- sample(test.input, shapiro.downsample.amount) %>% shapiro.test()
# transformed area
test.input <- halo.outputs.table[,normalized.area.variable]
shapiro.test.transformed <- sample(test.input, shapiro.downsample.amount) %>% shapiro.test()