#!/usr/bin/env R
#
# Get marker genes for cell types using mean ratio expression (log counts). 
# Uses DeconvoBuddies::get_mean_ratio2() function.
#

# devtools::install_github("https://github.com/LieberInstitute/DeconvoBuddies")
libv <- c("DeconvoBuddies", "SingleCellExperiment", "SummarizedExperiment")
sapply(libv, library, character.only = T)

#-------
# params
#-------
celltype.varname <- "cellType_broad_hc"

#----------
# set paths
#----------
proj.dpath <- "deconvo_method-paper"

# path to full singlecellexperiment
sce.fpath <- file.path(here(), 
                       "DLPFC_snRNAseq/processed-data/sce",
                       "sce_DLPFC.Rdata")

# save filepath
save.fnstem <- "k2"
sef.save.fname <- paste0("sef_mr-markers-",save.fnstem,
                         "-from-sce_dlpfc-ro1.rda")
markers.save.fname <- paste0("mr-markers-output_",save.fnstem,
                             "-ctbroadhc_from-sce_dlpfc-ro1.rda")
# make output filepaths
sef.fpath <- file.path(proj.dpath, "outputs/05_marker-gene-annotations/",
                       save.fname)
markers.fpath <- file.path(proj.dpath, "outputs/05_marker-gene-annotations/",
                           markers.save.fname)

#-----
# load
#-----
sce <- get(load(sce.fpath)) # get full singlecellexperiment

#---------------------
# summarize cell types
#---------------------
# variable name for cell types
table(sce[[celltype.varname]])
# Astro EndoMural     Micro     Oligo       OPC     Excit     Inhib
# 3979      2157      1601     32051      1940     24809     11067

#------------
# get markers
#------------
# get marker genes
markers <- get_mean_ratio2(sce, cellType_col = celltype.varname,
                           assay_name = "logcounts", add_symbol = TRUE)
# save results
save(markers, file = markers.fpath)

#------------------------------------------
# get sef -- subset top 100 markers by type
#------------------------------------------
table(markers$cellType.target) # num. markers returned by cell type
# Astro EndoMural     Micro     Oligo       OPC     Excit     Inhib
#   595       719       473       288      1333      4639      3676

# get top 100 markers by gene
ctv <- unique(markers$cellType.target)
ma <- as.data.frame(markers, stringsAsFactors = F)
ma[,2] <- as.character(ma[,2])
markerv <- unique(unlist(lapply(ctv, function(ki){
  mi <- as.data.frame(ma[ma$cellType.target == ki,])
  mi <- mi[order(mi$rank_ratio),]; mi[seq(100),1]
})))

# get counts for unique markers
scef <- sce[markerv,]
sef <- SummarizedExperiment(assays = list(counts = as.matrix(counts(scef))),
                            colData = colData(scef), 
                            rowData = rowData(scef))
# save sef
save(sef, file = sef.fpath)
