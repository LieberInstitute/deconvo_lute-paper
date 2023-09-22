library(ggplot2)

# load script 01 env
load("./env/07_summary/01_proportions_script.RData")

# load external proportions
sn.prop.external <- read.csv("./data/07_summary/snRNA_cell_type_proportions.csv")
halo.prop.external <- read.csv("./data/07_summary/HALO_cell_type_proportions.csv")

df.wide$neuron <- unlist(lapply(df.wide$sample.id, function(sample.id){
  filt.sn <- sn.prop.external$cell_type %in% c("Excit", "Inhib")
  filt.sn <- filt.sn & sn.prop.external$Sample==sample.id
  sum(sn.prop.external[filt.sn,]$prop_sn)
}))

ggplot(df.wide, aes(x = df.wide$neuron, y = df.wide$sn.neuron)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0)

ggplot(df.wide, aes(x = df.wide$neuron, y = df.wide$sn.neuron, color = sample.id)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0)

ggplot(df.wide[df.wide$sn.neuron<=0.2,], aes(x = neuron, y = sn.neuron, color = sample.id)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0) + xlim(0,1) + ylim(0,1)
