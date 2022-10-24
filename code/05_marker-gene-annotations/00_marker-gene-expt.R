#!/usr/bin/env R
#
# Get marker genes for cell types using mean ratio expression (log counts). 
# Uses DeconvoBuddies::get_mean_ratio2() function.
#

# devtools::install_github("https://github.com/LieberInstitute/DeconvoBuddies")
library(DeconvoBuddies) # contains get_mean_ratio2() to get marker genes
library(SingleCellExperiment)
library(SummarizedExperiment)
library(here)

#-------
# params
#-------
celltype.varname <- "cellType_broad_hc"

#----------
# set paths
#----------
proj.dpath <- "deconvo_method-paper"

#-----
# load
#-----
# path to full singlecellexperiment
sce.fpath <- file.path(here(), "DLPFC_snRNAseq/processed-data/sce/sce_DLPFC.Rdata")
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
markers <- get_mean_ratio2(sce, 
                           cellType_col = celltype.varname,
                           assay_name = "logcounts",
                           add_symbol = TRUE)
# save results
save(markers.withdrop, file = markers.withdrop.fpath)

#------------------------
# get markers --- no drop
#------------------------
# filter drop category
exclude.drop <- !sce[[celltype.varname]]=="drop"
scef <- sce[,exclude.drop]

# get marker genes
markers.nodrop <- get_mean_ratio2(scef, cellType_col = celltype.varname, 
                                  assay_name = "logcounts", add_symbol = TRUE)

# save results
save(markers.nodrop, file = markers.nodrop.fpath)

#----------------------
# marker gene summaries
#----------------------
lv <- list(markers.withdrop, markers.nodrop)
# number of unique marker genes
length(unique(markers.withdrop$gene)) # 8845
length(unique(markers.nodrop$gene)) # 5789

# top repeated marker genes
lapply(lv, function(ii){
  dt <- as.data.frame(table(ii$gene))
  dt <- dt[rev(order(dt[,2])),]
  print(head(dt))})
#        Var1 Freq
# 8757 ZNF638    8
# 8491 ZBTB20    8
# 8416   WWOX    8
# 7986   TTC3    8
# 7814 TNRC6B    8
# 7813 TNRC6A    8
#        Var1 Freq
# 5746 ZNF638    7
# 5616 ZBTB20    7
# 5562   WWOX    7
# 5274   TTC3    7
# 5163 TNRC6B    7
# 5162 TNRC6A    7