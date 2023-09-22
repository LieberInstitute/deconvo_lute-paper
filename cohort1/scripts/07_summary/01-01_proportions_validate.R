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
# check mae sce -- keep only overlapping cell type labels
# result: higher neuron proportions for sm
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

#--------------
# check mae sce -- remove ambiguous only
# result: closer than before, exact agreement for some sample IDs
#--------------
dim(sn)
celltypes.keep <- unique(sn.prop.external$cell_type)
filter.celltypes <- !sn[[celltype.variable.input]] %in% c("Ambiguous")
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


#--------------
# check mae sce -- use same script from : https://github.com/LieberInstitute/Human_DLPFC_Deconvolution/blob/02cd0b65950695c30dfd4426fedad7b730fc750c/code/03_HALO/08_explore_proportions.R#L74-L78
# result: very close proportions, comparable to shared proportions, use these ...
#--------------

library("tidyverse")
library("SingleCellExperiment")
library("ggrepel")
library("here")
library("sessioninfo")
library("here")
library("broom")

sce <- sn

load("./data/07_summary/cell_colors.RData")

halo_ct <- names(cell_type_colors_halo)

halo_ct_tb <- tibble(
  cell_type = c("Endo", "Astro", "Inhib", "Excit", "Micro", "Oligo"),
  # Other = rep(c("Other_Star", "Other_Circle"), each = 3),
  # marker = c("CLDN5", "GFAP", "GAD1", "SLC17A7", "TMEM119", "OLIG2"),
  Combo = rep(c("Circle", "Star"), each = 3)
)

sn_pd <- as.data.frame(colData(sce)) |>
  mutate(cell_type = factor(ifelse(gsub("Mural","",cellType_broad_hc) %in% halo_ct, 
                                   as.character(gsub("Mural","",cellType_broad_hc)),
                                   "Oligo"), ## OPC to Oligo
                            levels = halo_ct))

sn_pd |>
  dplyr::count(cellType_broad_hc, cell_type)

sn_ct <- sn_pd |> 
  select(Sample, cell_type) |> 
  left_join(halo_ct_tb) |> 
  bind_rows(sn_pd |> 
              select(Sample, cell_type) |> 
              left_join(tibble(
                cell_type = c("Endo", "Astro", "Inhib", "Excit", "Micro", "Oligo"),
                Combo = rep(c("Star", "Circle"), each = 3)
              ))|>
              mutate(cell_type = "Other")) |>
  as_tibble()

sn_ct |> dplyr::count(cell_type)

sn_ct_prop <- sn_ct |>
  group_by(Sample, cell_type, Combo) |>
  summarize(n_cell_sn = n()) |>
  group_by(Sample, Combo) |>
  mutate(prop_sn = n_cell_sn / sum(n_cell_sn)) 


df.wide$neuron.sn.sm.run <- unlist(lapply(df.wide$sample.id, function(sample.id){
  filt.sn <- sn_ct_prop$cell_type %in% c("Excit", "Inhib")
  filt.sn <- filt.sn & sn_ct_prop$Sample==sample.id
  sum(sn_ct_prop[filt.sn,]$prop_sn)
}))

ggplot(df.wide, aes(x = sn.neuron, y = neuron.sn.sm.run)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0)
