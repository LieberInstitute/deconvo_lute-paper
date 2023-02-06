#!/usr/bin/env R

#
#
#
#

libv <- c("scuttle", "dplyr", "limma", "ggplot2", "ggforce", "gridExtra",
          "glmGamPoi", "sva", "DeconvoBuddies", "SingleCellExperiment", "limma",
          "SummarizedExperiment")
sapply(libv, library, character.only = TRUE)

#----------
# load data
#----------
# get save dpath
code.dname <- "09_manuscript"
proj.dname <- "deconvo_method-paper"
save.dpath <- file.path(proj.dname, "outputs", code.dname)

# sce rda path
sce.fname <- "sce_DLPFC.Rdata"
sce.dpath <- "DLPFC_snRNAseq/processed-data/sce"
sce.fpath <- file.path(sce.dpath, sce.fname)

# get sce
# sce.fname <- "sce-mrb_dlpfc.rda"
# sce.fpath <- file.path(save.dpath, sce.fname)
# sce <- get(load(sce.fpath))

#-----------------------------------
# assign marker labels at variable k
#-----------------------------------
celltypevar <- "cellType_broad_hc"
table(sce[[celltypevar]])
# Astro EndoMural     Micro     Oligo       OPC     Excit     Inhib Ambiguous
# 3979      2157      1601     10894      1940     24809     11067     21157

# define marker categories
sce[["k2"]] <- ifelse(grepl("^Excit.*|^Inhib.*", sce[[celltypevar]]), "neuron", "other")
sce[["k3"]] <- ifelse(grepl("^Excit.*", sce[[celltypevar]]), "Excit", 
                      ifelse(grepl("^Inhib.*", sce[[celltypevar]]), "Inhib", "other"))
sce[["k4"]] <- ifelse(grepl("^Excit.*", sce[[celltypevar]]), "Excit", 
                      ifelse(grepl("^Inhib.*", sce[[celltypevar]]), "Inhib", 
                             ifelse(grepl("^Oligo$", sce[[celltypevar]]), "Oligo", "other")))

#----------------------------------
# global params for plots, analyses
#----------------------------------
# coldata
batchvar <- "donor"
marker.typev <- c("k2", "k3", "k4") # for iterations

# assays for adjustments
assayname.unadj <- "counts"
assayname.marker <- "logcounts"
assayname.adj <- "counts_adj"

# number of top markers per type
nmarker <- 20

# dispersion mean-var plots
smooth.linecol <- "cornflowerblue"
hl.color <- "red"
ref.linecol <- "black"

# dispersion from fits
num.genes.bg <- 1000

#------------------------
# run adjustment workflow
#------------------------
# format assay
assays(sce)[[assayname.unadj]] <- as.matrix(assays(sce)[[assayname.unadj]])

# run the workflow
for(markeri in marker.typev){
  message("working on marker: ", markeri, "...")
  message("doing combat adj...")
  mexpr <- assays(sce)[[assayname.unadj]]
  cnv <- colnames(sce)
  pheno <- data.frame(donor = sce[["donor"]], type = sce[[markeri]])
  mod <- model.matrix(~type, data = pheno)
  mi.adj <- ComBat(dat = mexpr, batch = pheno$donor, mod = mod)
  message("converting negative values...")
  mi.adj[mi.adj < 0] <- 0 # convert negative values
  assays(sce)[[assayname.adj]] <- mi.adj
  
  message("performing downsampling...")
  # downsample -- by donor, within celltypes
  utypev <- unique(sce[[markeri]])
  sce <- do.call(cbind, lapply(utypev, function(ti){
    message("downsampling for type ", ti, "...")
    # filter sce
    sce.filt <- sce[[markeri]]==ti
    scef <- sce[,sce.filt]
    # get filtered data
    batchv <- scef[[batchvar]]
    mexpr <- assays(scef)[[assayname.adj]]
    # downsample
    mexpr.ds <- downsampleBatches(mexpr, batch = batchv)
    assays(scef)[[assayname.adj]] <- mexpr.ds 
    scef
  }))
  
  message("Getting top markers..")
  sce <- logNormCounts(sce, assay.type = assayname.adj)
  mr <- get_mean_ratio2(sce, assay_name = "logcounts", cellType_col = markeri)
  # get top N markers from results
  typev <- unique(mr$cellType.target)
  mrtop <- do.call(rbind, lapply(typev, function(typei){
    mr %>% filter(cellType.target == typei) %>% 
      arrange(rank_ratio) %>% top_n(n = nmarker)
  }))
  # store in sce
  metadata(sce)[["markers"]] <- list(all = mr, top = mrtop)
  
  # filter assays
  sce.filt <- names(assays(sce)) %in% c(assayname.unadj, assayname.adj)
  assays(sce) <- assays(sce)[sce.filt]
  
  # append metadata
  metadata(sce)[["adj_method"]] <- "combat"
  metadata(sce)[["ds_method"]] <- list(package = "scuttle",
                                       funct = "downsampleBatches",
                                       typevar = markeri, groupvar = "donor")
  
  # save
  fname <- paste0("sce_marker-adj-",markeri,"_mrb-dlpfc.rda")
  save(sce, file = file.path(save.dpath, fname))
  message("finished with marker ", markeri)
}

# returns:
# working on marker: k2...
# doing combat adj...
# Found 7406 genes with uniform expression within a single batch (all zeros); these will not be adjusted for batch.
# ...
# working on marker: k3...
# doing combat adj...
# Found 7406 genes with uniform expression within a single batch (all zeros); these will not be adjusted for batch.
# ...
