#!/usr/bin/env R

# Author: Sean Maden
#

source("deconvo_method-paper/code/_halo-outputs-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
halo.outputs.table <- get(load(halo.output.path))
halo.outputs.table <- halo.outputs.table %>% as.data.frame()
anova.result.list <- get(load(anova.results.list.path))

# make residual plots
# get residuals
residuals1 <- residuals(anova.model1)
residuals2 <- residuals(anova.model2)
residuals3 <- residuals(anova.model3)
# do quantile-quantile plots
# do quantile-quantile plots
qqnorm(residuals1)
qqnorm(residuals2)
qqnorm(residuals3)
# do quantile-quantile lines
qqline(residuals1)
qqline(residuals2)
qqline(residuals3)