#!/usr/bin/env R

# Author: Sean Maden
#
#

libv <- c("SingleCellExperiment", "SummarizedExperiment")
sapply(libv, library, character.only = T)

#----------
# load data
#----------
# load bulk data
rse.filename <- "rse_gene.Rdata"
load.path <- file.path("Human_DLPFC_Deconvolution", "processed-data",
                       "01_SPEAQeasy", "round2_v40_2022-07-06", "rse")
rse <- get(load(file.path(load.path, rse.filename)))

# get save directory path
save.path <- file.path("deconvo_method-paper", "outputs", 
                       "11_bulk-expression-analysis")

# load marker data
sce.markers.path <- file.path("deconvo_method-paper", "outputs", "09_manuscript", 
                              "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda")
lsce <- get(load(sce.markers.path))
sce <- lsce[["k2"]]
rm(lsce)

# load helper functions
group.jitter.filename <- "group-jitter_helper.rda"
get.summary.filename <- "get-summary_helper.rda"
group_jitter <- get(load(file.path(save.path, group.jitter.filename)))
get_summary_list <- get(load(file.path(save.path, get.summary.filename)))

#------------------------
# manage helper functions
#------------------------

#--------------------------------
# get marker gene bulk expression
#--------------------------------
# subset rse
bulk.gene.names <- rowData(rse)$Symbol
marker.genes.vector <- rownames(sce)
overlapping.markers <- intersect(bulk.gene.names, marker.genes.vector)
message("Found ", length(overlapping.markers), " overlapping markers.")
filter <- which(rowData(rse)$Symbol %in% overlapping.markers)
rsef <- rse[filter,]
dim(rsef)

# save bulk marker expr
rsef.filename <- "rsef_k2-marker-expr_ro1-dlpfc.rda"
save.path <- file.path(save.path, rsef.filename)
save(rsef, file = save.path)

#--------------------------------
# qc at markers versus background
#--------------------------------
# set up data summaries
# params
assay.name <- "counts"
batch.variable <- "batch.id"
condition.variable <- "expt_condition"
# rse data
cd <- colData(rse)
counts.bg <- assays(rse)[["counts"]]
counts.marker <- assays(rsef)[["counts"]]
# get new cd variables
cd[,batch.variable] <- paste0(cd$BrNum,"_",cd$location)
cd[,condition.variable] <- paste0(cd$library_prep,"_",cd$library_type)
# vector of group variables to summarize
variable.vector <- c(batch.variable, "library_prep", "library_type", 
                     condition.variable)
# define summary types
type.vector <- c("total.counts", "dispersion", "zero.count", "mean", "variance")

for(type in type.vector){
  message("Working on summary type ", type, "...")
  
  message("Working on background expression...")
  
  message("Working on marker expression...")
  
  message("Getting new expression data...")
  
  get_summary_list(type = type, plot.fname = plot.fname, 
                   variable.vector = variable.vector, cd = cd, 
                   counts = counts, save.path = save.path)
  message("Finished with summary type ", type, ".")
  
}
message("Done with all summary  types.")








