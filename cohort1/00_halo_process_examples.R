

load("./data/09_quality/halo_all.Rdata")
load("./env/03_shuffle/00_fig3ab_script.RData")

halo.table = halo_all
data.type = "halo_all"
sample.id.colname = "Sample"
combo.label.varname = "Combo"
combo.label.vector = c("Circle", "Star")

halo.processed.table