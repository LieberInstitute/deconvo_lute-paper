#!/usr/bin/env R

# Author: Sean Maden
#
# Summarize results of within matched samples deconvolution experiments.

source("deconvo_method-paper/code/14_deconvolution-framework-trials/00_parameters.R")
sapply(libv, library, character.only = T)
results.table <- get(load(within.samples.results.table.path))
results.table$method <- gsub("Param", "", results.table$method.string)

# compare total cells
results.table$total.cells <- results.table$glial.count.true+results.table$neuron.count.true
summary(results.table$total.cells)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 9461   11189   27524   27970   41740   53710
median(results.table$total.cells) # 27524
filter <- results.table$absolute.error.neuron < 0.2
median(results.table[filter,]$total.cells) # 20019
filter <- results.table$absolute.error.neuron >= 0.2
median(results.table[filter,]$total.cells) # 27524

# compare total expression
summary(results.table$y.total.expression)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 15773  207003  318027  344890  401559 1272108
filter <- results.table$absolute.error.neuron < 0.1
summary(results.table[filter,]$y.total.expression) 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 166636  184824  279416  519201  815898 1272108
filter <- results.table$absolute.error.neuron >= 0.1
summary(results.table[filter,]$y.total.expression)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 15773  207003  330492  309382  401559  587330

# expression bulk by cell count image
ggplot(results.table, aes(x = total.cells, y = y.total.expression,
                          size = absolute.error.neuron)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0)
# exaggerated errors
results.table$abs.error.sq <- (results.table$absolute.error.neuron*10)^2
ggplot(results.table, aes(x = total.cells, y = y.total.expression,
                          size = abs.error.sq)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0)

# correlate, total cells by total bulk expression
cor.test(results.table$total.cells, 
         results.table$y.total.expression)

