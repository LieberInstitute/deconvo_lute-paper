#!/usr/bin/env R

# Author: Sean Maden
#
# Plot results of within matched samples deconvolution experiments.

source("deconvo_method-paper/code/14_deconvolution-framework-trials/00_parameters.R")
sapply(libv, library, character.only = T)
results.table <- get(load(within.samples.results.table.path))

# plot cell proportions
new.plot <- ggplot(results.table, 
                   aes(x = neuron.proportion.true, y = neuron_proportion,
                       color = method)) +
  geom_point(alpha = 0.4) + theme_bw() +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  geom_vline(xintercept = 0.5, color = "black", alpha = 0.4) +
  xlab("True") + ylab("Predicted") + ggtitle("Neuron cell proportion") +
  xlim(0, 1) + ylim(0, 1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#
jpeg(scatterplot.proportions.bysample.colmethod.path, width = 7, height = 5, 
     units = "in", res = 400)
new.plot + facet_wrap(~sample.id)
dev.off()
#
jpeg(scatterplot.proportions.bylibprep.colmethod.path, width = 6, height = 2.5, 
     units = "in", res = 400)
new.plot + facet_wrap(~library.prep)
dev.off()
#
jpeg(scatterplot.proportions.bylibtype.colmethod.path, width = 5, height = 2.5,
     units = "in", res = 400)
new.plot + facet_wrap(~library.type)
dev.off()
#
jpeg(scatterplot.proportions.byexptgroup.colmethod.path, width = 6, height = 3.5,
     units = "in", res = 400)
new.plot + facet_wrap(~bulk.id)
dev.off()
#
jpeg(scatterplot.proportions.bymethod.colmethod.path, width = 5.2, height = 2,
     units = "in", res = 400)
new.plot + facet_wrap(~method)
dev.off()

# plot error
new.plot <- ggplot(results.table, aes(x = method, y = absolute.error.neuron)) + theme_bw() + 
  geom_jitter(alpha = 0.5) + geom_boxplot(alpha = 0, color = "cyan") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0,1) +
  ylab("Absolute error") + ggtitle("Cell type: neuron") + 
  xlab("Deconvolution algorithm")
#
jpeg(jitterbox.abserror.bysample.xmethod.path, width = 6, height = 4, 
     units = "in", res = 400)
new.plot + facet_wrap(~sample.id)
dev.off()
#
jpeg(jitterbox.abserror.bylibprep.xmethod.path, width = 6, height = 4, 
     units = "in", res = 400)
new.plot + facet_wrap(~library.prep)
dev.off()
#
jpeg(jitterbox.abserror.bylibtype.xmethod.path, width = 6, height = 4, 
     units = "in", res = 400)
new.plot + facet_wrap(~library.type)
dev.off()
#
jpeg(jitterbox.abserror.byexptgroup.xmethod.path, width = 6, height = 4, 
     units = "in", res = 400)
new.plot + facet_wrap(~bulk.id)
dev.off()
#

# plot rmse
new.plot <- ggplot(results.table, aes(x = method, y = rmse)) + theme_bw() + 
  geom_jitter(alpha = 0.2) + geom_boxplot(alpha = 0, color = "cyan") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Deconvolution algorithm") + ylab("RMSE") + ggtitle("Cell type: neuron")
#
jpeg(jitterbox.rmse.bysample.xmethod.path, width = 6, height = 4, 
     units = "in", res = 400)
new.plot + facet_wrap(~sample.id)
dev.off()
#
jpeg(jitterbox.rmse.bylibprep.xmethod.path, width = 5, height = 3, 
     units = "in", res = 400)
new.plot + facet_wrap(~library.prep)
dev.off()
#
jpeg(jitterbox.rmse.bylibtype.xmethod.path, width = 5, height = 3, 
     units = "in", res = 400)
new.plot + facet_wrap(~library.type)
dev.off()
#
jpeg(jitterbox.rmse.byexptgroup.xmethod.path, width = 5, height = 4, 
     units = "in", res = 400)
new.plot + facet_wrap(~bulk.id)
dev.off()
