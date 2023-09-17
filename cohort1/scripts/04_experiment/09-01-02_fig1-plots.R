# fig 1a plots
experiment.name <- "experiment_within-slide-matched_dlpfc-ro1-train.rda"
experiment.path <- here("deconvo_method-paper", "outputs", experiment.name)
experiment <- get(load(experiment.path))
# plot
plot.data <- experiment$results.table
plot.data <- plot.data[plot.data$cell.size.adjustment.type=="null",]
jpeg("fig1a.jpg", width = 10, height = 10, units = "in", res = 400)
ggplot(plot.data, aes(x = neuron.true.proportion, y = neuron.predicted.proportion)) +
  geom_point() + geom_abline(slope = 1, intercept = 0) + theme_bw() +
  geom_hline(yintercept = 0.5, alpha = 0.5) + 
  geom_vline(xintercept = 0.5, alpha = 0.5) +
  xlab("True proportion") + ylab("Predicted proportion") +
  xlim(0, 1) + ylim(0, 1)
dev.off()