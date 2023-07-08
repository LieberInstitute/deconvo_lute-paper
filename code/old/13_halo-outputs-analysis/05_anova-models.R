#!/usr/bin/env R

# Author: Sean Maden
#
# Perform analyses of variance (ANOVAs).
#

source("deconvo_method-paper/code/13_halo-outputs-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
halo.outputs.table <- get(load(output.updated.path))
halo.outputs.table <- halo.outputs.table %>% as.data.frame()

# set dependent variable
halo.outputs.table$anova_depenent_variable <- halo.outputs.table[,normalized.area.variable]

# set models
basic.model.string <- "anova_depenent_variable ~ cell_type + Slide + SAMPLE_ID"
complex.model.string <- "anova_depenent_variable ~ cell_type + Slide + Combo + Position + BrNum"

# linear models
# basic
basic.linear.model.string <- paste0("lm(", basic.model.string, ", data = halo.outputs.table)")
basic.linear.model <- eval(parse(text = basic.linear.model.string))
# complex
complex.linear.model.string <- paste0("lm(", complex.model.string, ", data = halo.outputs.table)")
complex.linear.model <- eval(parse(text = complex.linear.model.string))

# anova tests
# basic
anova.string <- paste0("aov(", basic.model.string, ", data = halo.outputs.table)")
basic.anova.test <- eval(parse(text = anova.string))
# complex
anova.string <- paste0("aov(", complex.model.string, ", data = halo.outputs.table)")
complex.anova.test <- eval(parse(text = anova.string))

# save anova results list
model.results.list <- list(basic = list(linear.model = basic.linear.model, 
                                        anova = basic.anova.test),
                           complex = list(linear.model = basic.linear.model, 
                                          anova = complex.anova.test))
save(model.results.list, file = model.results.list.path)
