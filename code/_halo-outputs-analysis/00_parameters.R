libv <- c("nlme", "ggplot2", "gridExtra", "dplyr")
sapply(libv, library, character.only = TRUE)

# set the save directory
save.path <- here("deconvo_method-paper", "outputs", "09_manuscript")
# set the halo data path
halo.output.file.name <- "halo_all.Rdata"
halo.output.path <- here("Human_DLPFC_Deconvolution", "processed-data", 
                              "03_HALO", halo.output.file.name)

# cell labels
labels <- c("Endo" = "CLDN5", "Astro" = "GFAP", "Inhib" = "GAD1", 
            "Excit" = "SLC17A7", "Micro" = "TMEM119", "Oligo" = "OLIG2")

# 01
sample.id.label <- "Sample"
cell.area.variable <- "Nucleus_Area"
cell.area.log.variable <- "log10_nucleus_area"
gene.marker.label <- "AKT3_Copies"
output.updated.filename <- "halo_updated_path.Rdata"
output.updated.path <- here(save.path, output.updated.filename)
marker.quantile.variable <- "akt3.copies.quantile.scale"
halo.quantiles.jpg.file.name <- ""

# 02

# 99 quantile scale summaries
