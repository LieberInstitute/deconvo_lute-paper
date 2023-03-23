#!/usr/bin/env R

# Author: Sean Maden
#

source("deconvo_method-paper/code/_halo-outputs-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
halo.outputs.table <- get(load(halo.output.path))
halo.outputs.table <- halo.outputs.table %>% as.data.frame()

# get anova models
model.string <- paste0(anova.dependent.variable, " ~ cell_type + Slide + SAMPLE_ID")
model1 <- paste0("lm(", model.string, ", data = halo.outputs.table)") %>% 
  parse() %>% eval()
model2 <- paste0("gls(", model.string, ", data = halo.outputs.table, method = 'REML')") %>% 
  parse() %>% eval()
model3 <- paste0("glm(", model.string, ", data = halo.outputs.table)") %>% 
  parse() %>% eval()

# anova analysis
dependent.variable <- cell.area.log.variable
# basic models
model.string <- paste0(dependent.variable, " ~ cell_type + Slide + SAMPLE_ID")
anova.string <- paste0("aov(", model.string, ", data = halo.outputs.table)")
anova.model1 <- anova.string %>% parse() %>% eval()
# complex model
model.string <- paste0(dependent.variable, " ~ cell_type + Slide + Combo + Position + BrNum")
anova.string <- paste0("aov(", model.string, ", data = halo.outputs.table)")
anova.model2 <- anova.string %>% parse() %>% eval()
# akt3
dependent.variable <- gene.marker.label
model.string <- paste0(dependent.variable, " ~ cell_type + Slide + Combo + Position + BrNum")
anova.test.string <- paste0("aov(", model.string, ", data = halo.outputs.table)")
result1 <- anova.test.string %>% parse() %>% eval() # eval(parse(text = anova.test.string))
# complex model
model.string <- paste0(dependent.variable, " ~ cell_type + Slide + Combo + Position + BrNum")
anova.string <- paste0("aov(", model.string, ", data = halo.outputs.table)")
result2 <- anova.string %>% parse() %>% eval()
