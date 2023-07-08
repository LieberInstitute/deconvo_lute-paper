#!/usr/bin/env R

#
# Filter the bulk RNA-seq data on gene types and k2 markers.
#

# load data
libv <- c("here", "SummarizedExperiment", "SingleCellExperiment", "dplyr", "data.table", "DESeq2", "ggplot2", "ggcorrplot", "glmGamPoi", "gridExtra")
sapply(libv, library, character.only = T)

# paths
save.path <- here("deconvo_method-paper", "outputs", "01_prepare-datasets")
rse.bulk.filename <- "rse_gene.Rdata"
rse.bulk.filepath <- here("Human_DLPFC_Deconvolution", "processed-data", "01_SPEAQeasy", "round2_v40_2022-07-06", "rse", rse.bulk.filename)
rse.bulk.filename.new <- "rse-bulk-resave.rda"
rse.bulk.path.new <- here(save.path, rse.bulk.filename.new)
rse.k2markers.filename <- "rse_k2-marker-expression_ro1-dlpfc.rda"
rse.k2markers.filepath <- here(save.path, rse.k2markers.filename)
rse.gene.filter.filename <- "rse-gene-filter.rda"
rse.gene.filter.filepath <- here(save.path, rse.gene.filter.filename)
sce.markers.filename <- "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda"
sce.markers.list.path <- here(save.path, sce.markers.filename)

# variable names
experiment.condition1 <- "library_prep"
experiment.condition2 <- "library_type"
condition.variable <- "expt_condition"
donor.variable <- "BrNum"
location.variable <- "location"
batch.variable <- "batch.id"
# gene type variables
gene.types.protein <- c("protein_coding")
gene.types.nonpolya <- c("lncRNA", "Mt_rRNA", "rRNA", "Mt_tRNA")
gene.types.include <- c(gene.types.protein, gene.types.nonpolya)

# load
rse <- get(load(file.path(rse.bulk.filepath)))

# view gene types
rd <- rowData(rse)
table(rd$gene_type)
cd <- colData(rse)
# add batch id
cd[,batch.variable] <- paste0(cd[,donor.variable], "_", cd[,location.variable])
# add experiment groups
cond1 <- cd[,experiment.condition1]
cond2 <- cd[,experiment.condition2]
cd[,condition.variable] <- paste0(cond1, "_", cond2)
# reassign coldata
colData(rse) <- cd

# filter gene types
gene.types.vector <- rd$gene_type
gene.type.filter <- which(gene.types.vector %in% gene.types.include) 
rse.filter <- rse[gene.type.filter,]

# save
# save(rse, file = rse.bulk.path.new)
save(rse.filter, file = rse.gene.filter.filepath)

# subset k2 markers
# load k2 sce data
sce <- get(load(sce.markers.list.path))[["k2"]]
# get bulk marker expression
marker.genes.vector <- rownames(sce)
# subset rse on markers
bulk.gene.names <- rowData(rse.filter)$Symbol
gene.marker.intersect <- intersect(bulk.gene.names, marker.genes.vector)
message("Found ", length(gene.marker.intersect), " overlapping markers.")
rse.filter.markers <- which(bulk.gene.names %in% gene.marker.intersect)
rse.filter.markers <- rse.filter[rse.filter.markers,]
# gene symbols as rownames
rownames(rse.filter.markers) <- rowData(rse.filter.markers)$Symbol
dim(rse.filter.markers)
# save
save(rse.filter.markers, file = rse.k2markers.filepath)