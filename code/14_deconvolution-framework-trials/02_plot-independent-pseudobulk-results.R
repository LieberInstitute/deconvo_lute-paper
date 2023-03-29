#!/usr/bin/env R

# Author: Sean Maden
#
# Plot independent pseudobulk results.

source("deconvo_method-paper/code/14_deconvolution-framework-trials/00_parameters.R")
sapply(libv, library, character.only = T)
results.table <- get(load(independent.pb.results.table.path))

# plot true vs. pred proportions, neuron
new.plot <- ggplot(results.table, aes(x = neuron.true.prop, 
                                      y = neuron.prop.pred, 
                                      color = method)) +
  geom_point(alpha = 0.6) + geom_abline(slope = 1, intercept = 0, color = "black") +
  geom_vline(xintercept = 0.5, alpha = 0.5) + xlab("True") + ylab("Predicted") +
  ggtitle("Cell type: neuron") + xlim(0, 1) + ylim(0, 1) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

new.plot + facet_wrap(~sample.id)
jpeg(pb.scatterplot.proportions.bysample.colmethod.path, width = 5.5, 
     height = 2.3, units = "in", res = 400)
new.plot + facet_wrap(~sample.id)
dev.off()

# plot absolute errors, neuron
results.table$neuron.abs.error <- abs(results.table$neuron.error)
new.plot <- ggplot(results.table, 
                   aes(x = method, y = neuron.abs.error, fill = method)) +
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  xlab("Deconvolution algorithm") +
  ylab("Absolute error") +
  ggtitle("Cell type: neuron")
jpeg(pb.barplot.abserror.bysample.colmethod.path, width = 5, height = 4, 
     units = "in", res = 400)
new.plot + facet_wrap(~sample.id)
dev.off()

# plot rmse, k2
new.plot <- ggplot(results.table, aes(x = method, y = rmse.k2, fill = method)) +
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  xlab("Deconvolution algorithm") +
  ylab("RMSE") +
  ggtitle("Cell type: neuron")
jpeg(pb.barplot.rmse.bysample.colmethod.path, width = 5, height = 4, 
     units = "in", res = 400)
new.plot + facet_wrap(~sample.id)
dev.off()
