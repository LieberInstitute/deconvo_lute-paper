#!/usr/bin/env R

# Author: Sean Maden
#

source("deconvo_method-paper/code/_halo-outputs-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
halo.outputs.table <- get(load(halo.output.path))
halo.outputs.table <- halo.outputs.table %>% as.data.frame()

# explained variance analysis
avi <- results1
dfp <- summary(avi)[[1]] %>% as.data.frame()
dfp <- data.frame(variable = gsub(" ", "", rownames(dfp)), ssq = dfp[,2])
dfp$perc.var <- 100*dfp$ssq/sum(dfp$ssq)
dfp$label <- paste0(round(dfp$perc.var, 1), "%")
level.order <- dfp$perc.var %>% order() %>% rev()
dfp$variable <- factor(dfp$variable, levels = dfp$variable[level.order])

# bar plot with residuals
ggbar <- ggplot(dfp, aes(x = variable, y = perc.var, fill = variable)) +
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Total percent variance") + ggtitle(title.str) +
  geom_text(aes(x = variable, y = perc.var, label = label)) +
  theme(legend.position = "none")
jpeg(barplot.residuals.jpg.path, width = 3.5, height = 2.5, units = "in", res = 400); ggbar; dev.off()

# bar plot with residuals, log10 scale
ggbar <- ggplot(dfp, aes(x = variable, y = perc.var, fill = variable)) +
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  ylab("Percent variance\n(log10-scaled)") + 
  scale_y_log10() + ggtitle(title.str) +
  geom_text(aes(x = variable, y = perc.var, label = label))
jpeg(barplot.residuals.log10.jpg.path, width = 3.5, height = 2.5, units = "in", res = 400); ggbar; dev.off()

# bar plot without residuals
dfpe <- dfp[!dfp$variable=="Residuals",]
dfpe$perc.var <- 100*dfpe$ssq/sum(dfpe$ssq)
dfpe$label <- paste0(round(dfpe$perc.var, 1), "%")
ggbar <- ggplot(dfpe, aes(x = variable, y = perc.var, fill = variable)) +
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  ylab("Total explained variance") + ggtitle(title.str) +
  geom_text(aes(x = variable, y = perc.var, label = label))
jpeg(barplot.noresiduals.jpg.path, width = 3.5, height = 2.5, units = "in", res = 400); ggbar; dev.off()
