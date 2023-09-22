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


#--------------
# check mae sce
#--------------
mae <- get(load("./outputs/01_mae/mae_allsamples.rda"))
sn <- mae[["snrnaseq.k2.all"]]

sn.prop.external <- read.csv("./data/07_summary/snRNA_cell_type_proportions.csv")
unique.cell.types <- unique(sn.prop.external$cell_type)

celltype.variable.input <- "cellType_broad_hc"
dim(sn)
celltypes.keep <- unique(sn.prop.external$cell_type)
filter.celltypes <- sn[[celltype.variable.input]] %in% celltypes.keep
sn.new <- sn[,filter.celltypes]
dim(sn.new)
celltypevar <- "k2"
sn.new[[celltypevar]] <- ifelse(
  grepl("^Excit.*|^Inhib.*", sn.new[[celltype.variable.input]]), "neuron", "glial")

df.prop <- se_cell_prop(sn.new, "k2", "Sample", label = "revised")
colnames(df.prop) <- paste0(colnames(df.prop), ".sn.new")

df.prop$neuron.external <- unlist(lapply(df.prop$sample.id, function(sample.id){
  filt.sn <- sn.prop.external$cell_type %in% c("Excit", "Inhib")
  filt.sn <- filt.sn & sn.prop.external$Sample==sample.id
  sum(sn.prop.external[filt.sn,]$prop_sn)
}))

ggplot(df.prop, aes(x = neuron.external, y = neuron.sn.new)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0)






