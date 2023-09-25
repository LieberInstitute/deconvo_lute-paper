
# load sce
mae <- get(load("./outputs/01_mae/mae_allsamples.rda"))
sn <- mae[["snrnaseq.k2.all"]]

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

load("./data/cell_colors.RData")

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

