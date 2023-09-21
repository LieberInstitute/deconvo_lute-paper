#!/usr/bin/env R

#
# Get cell type proportions
#

libv <- c("dplyr", "ggplot2")
sapply(libv, library, character.only = TRUE)

source("./scripts/07_summary/00_prop.R")

# load
mae <- get(load("./outputs/01_mae/mae_allsamples.rda"))
sn <- mae[["snrnaseq.k2.all"]]
sample.id.variable <- unique(sn[["Sample"]])
mae <- mae[,colData(mae)$sample.id %in% sample.id.variable,]

#----------------
# get proportions
#----------------
list.dfp <- list_dfp_wide_tall(mae)

# save env
save.image("./env/07_summary/01_proportions_script.RData")

#------
# plots
#------
dfp.wide <- list.dfp$dfp.wide
dfp.tall <- list.dfp$dfp.tall

# wide plots
ggplot(dfp.wide, aes(x = sn.neuron, y = rn.neuron)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0)

ggplot(dfp.wide, aes(x = sn.neuron, y = rn.neuron)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0) + facet_wrap(~sn.sample.id)

ggplot(dfp.wide, aes(x = sn.neuron, y = rn.neuron)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0) + facet_wrap(~br.region)

ggplot(dfdfp.wide, aes(x = sn.neuron, y = rn.neuron)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0) + facet_wrap(~subject.id)

# tall plots
ggplot(dfp.tall, aes(x = prop.type.label, y = neuron)) + 
  geom_violin(draw_quantiles = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(dfp.tall, aes(x = prop.type.label, y = neuron)) + 
  geom_violin(draw_quantiles = 0.5) + facet_wrap(~br.region) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(dfp.tall, aes(x = prop.type.label, y = neuron)) + 
  geom_violin(draw_quantiles = 0.5) + facet_wrap(~subject.id) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(dfp.tall, aes(x = prop.type.label, y = neuron)) + 
  geom_jitter() + geom_boxplot(alpha = 0, color = "cyan") + facet_wrap(~subject.id) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(dfp.tall, aes(x = prop.type.label, y = neuron)) + 
  geom_jitter() + geom_boxplot(alpha = 0, color = "cyan") + facet_wrap(~br.region) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

###