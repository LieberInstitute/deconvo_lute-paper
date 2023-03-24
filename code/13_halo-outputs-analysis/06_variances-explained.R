#!/usr/bin/env R

# Author: Sean Maden
#

source("deconvo_method-paper/code/13_halo-outputs-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
model.results.list <- get(load(model.results.list.path))

# basic model results
anova.results <- anova.results <- model.results.list$basic$anova
plot.data.frame <- summary(anova.results)[[1]] %>% as.data.frame()
plot.data.frame <- data.frame(variable = gsub(" ", "", rownames(plot.data.frame)), 
                              ssq = plot.data.frame[,2])
plot.data.frame$perc.var <- 100*plot.data.frame$ssq/sum(plot.data.frame$ssq)
plot.data.frame$label <- paste0(round(plot.data.frame$perc.var, 1), "%")
level.order <- plot.data.frame$perc.var %>% order() %>% rev()
plot.data.frame$variable <- factor(plot.data.frame$variable, 
                                   levels = plot.data.frame$variable[level.order])
# bar plot with residuals
ggplot.barplot <- ggplot(plot.data.frame, aes(x = variable, y = perc.var, fill = variable)) +
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Total percent variance") + ggtitle(plot.title.string.basic) +
  geom_text(aes(x = variable, y = perc.var, label = label)) +
  theme(legend.position = "none") + facet_zoom(ylim = c(0,3))
# save
jpeg(barplot.basic.variance.explained.path, 
     width = 8, height = 2.5, units = "in", res = 400)
ggplot.barplot
dev.off()

# complex model results
anova.results <- anova.results <- model.results.list$complex$anova
plot.data.frame <- summary(anova.results)[[1]] %>% as.data.frame()
plot.data.frame <- data.frame(variable = gsub(" ", "", rownames(plot.data.frame)), 
                              ssq = plot.data.frame[,2])
plot.data.frame$perc.var <- 100*plot.data.frame$ssq/sum(plot.data.frame$ssq)
plot.data.frame$label <- paste0(round(plot.data.frame$perc.var, 1), "%")
level.order <- plot.data.frame$perc.var %>% order() %>% rev()
plot.data.frame$variable <- factor(plot.data.frame$variable, 
                                   levels = plot.data.frame$variable[level.order])
# bar plot with residuals
ggplot.barplot <- ggplot(plot.data.frame, aes(x = variable, y = perc.var, fill = variable)) +
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Total percent variance") + ggtitle(plot.title.string.complex) +
  geom_text(aes(x = variable, y = perc.var, label = label)) +
  theme(legend.position = "none") + facet_zoom(ylim = c(0,3))
# save
jpeg(barplot.complex.variance.explained.path, 
     width = 8, height = 2.5, units = "in", res = 400)
ggplot.barplot
dev.off()
