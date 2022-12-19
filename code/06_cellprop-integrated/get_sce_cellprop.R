#!/usr/bin/env R

# Author: Sean Maden
#
# Get cell proportions from main SCE object

libv <- c("SingleCellExperiment", "SummarizedExperiment")
sapply(libv, library, character.only = T)

#----------
# load data
#----------
# load sce
sce.fname <- "sef_mr-markers-k7-from-sce_dlpfc-ro1.rda"
sce.fpath <- file.path("deconvo_methods-paper", "outputs", "05_marker-gene-annotations", sce.fname)
sce <- get(load(sce.fpath))

# load halo
halo.fpath <- file.path("Human_DLPFC_Deconvolution", "raw-data", "HALO", 
                        "Algorithm_Check_20220914", "Halo_CSV")
halo.lfv <- list.files(halo.fpath)
halo.lfv <- halo.lfv[grepl(".*\\.csv$", halo.lfv)]
halo.lfv.paths <- file.path(halo.fpath, halo.lfv)
lh <- lapply(halo.lfv.paths, read.csv); names(lh) <- halo.lfv
# load rnascope data
fpath <- file.path("Human_DLPFC_Deconvolution", "raw-data", "HALO",
                   "Algorithm_Check_20220914", "Halo_CSV")
fnv <- list.files(fpath, pattern = ".*.csv$", recursive=T) # get long paths
fnv <- fnv[!grepl("prelim", fnv)]; dirv <- c("CIRCLE", "STAR")
lcsv <- lapply(dirv, function(diri){
  fnvi <- fnv[grepl(paste0(".*",diri,".*"), fnv)]
  lcsvi <- lapply(fnvi, function(fni){read.csv(file.path(fpath, fni))})
  names(lcsvi) <- gsub(".*\\/|\\.csv$", "", fnvi)
  return(lcsvi)
})
names(lcsv) <- dirv

#-----------------
# helper functions
#-----------------
cellprop_from_halo <- function(csvi, samplab, region.str, expt = c("CIRCLE", "STAR"),
                               labexpt = list(CIRCLE = c("Endo", "Astro", "Inhib"),
                                              STAR = c("Excit", "Micro", "Oligo")),
                               labels = c("Endo" = "CLDN5", "Astro" = "GFAP",
                                          "Inhib" = "GAD1", "Excit" = "SLC17A7",
                                          "Micro" = "TMEM119", "Oligo" = "OLIG2")){
  # get labels for experiment
  labf <- labels[labexpt[[expt]]]
  # check cols
  cnvf <- which(grepl(paste(paste0("^", labf, "$"), collapse = "|"),
                      colnames(csvi)))
  dfi <- csvi[,cnvf]
  # relabel multi-mappers
  rsumv <- rowSums(dfi); which.mm <- which(rsumv > 1)
  dfi[which.mm,1] <- dfi[which.mm,2] <- dfi[which.mm,3] <- 0 
  # get proportions
  do.call(rbind, lapply(seq(length(labf)), function(ii){
    celli <- names(labf)[ii]; markeri <- labf[[ii]]
    num.celli <- length(which(dfi[,markeri]==1))
    data.frame(cell_type = celli, num_cells = num.celli, 
               prop_cells = num.celli/nrow(dfi), region = region.str, 
               sample_id = samplab, expt = expt)
  }))
}


#---------
# get prop
#---------
# prop -- sce
donor.varname.sce <- "BrNum"
ct.varname.sce <- "cellType_broad_hc"
dv <- unique(sce[[donor.varname.sce]])
dfp <- do.call(rbind, lapply(dv, function(di){
  message(di)
  scef <- sce[,sce[[donor.varname.sce]]==di]
  ctvi <- unique(scef[[ct.varname.sce]])
  dfi <- do.call(rbind, lapply(ctvi, function(ci){
    ncol(scef[,scef[[ct.varname.sce]]==ci])}))
  colnames(dfi) <- "cells"
  dfi <- as.data.frame(dfi)
  dfi$prop <- dfi$cells/sum(dfi$cells)
  dfi$cell.type <- ctvi
  dfi$donor <- di
  return(dfi)
}))

# prop --- halo
# rnascope slides
sampv.rn <- unique(sapply(names(lcsv[[1]]), function(ni){
  stri <- strsplit(ni, "_")[[1]][3]; gsub("[A-Z]", "", stri)}))

dfrn <- do.call(rbind, lapply(c("CIRCLE", "STAR"), function(expt){
  lcsve <- lcsv[[expt]]
  do.call(rbind, lapply(unique(dfsn$sample_id), function(sampi){
    lcsvi <- lcsve[grepl(paste0("_", sampi, "[A-Z]"), names(lcsve))]
    if(length(lcsvi) > 0){
      do.call(rbind, lapply(seq(length(lcsvi)), function(ii){
        # get region string
        namei <- names(lcsvi)[ii]
        region.str <- gsub("[0-9]", "", strsplit(namei, "_")[[1]][3])
        region.str <- ifelse(region.str == "M", "middle",
                             ifelse(region.str == "A", "anterior", "posterior"))
        # make new cell prop df
        cellprop_from_halo(lcsvi[[ii]], sampi, region.str, expt)
      }))
    }
  }))
}))


# prop --- bulk

#-------------
# analyze prop
#-------------
# correlations
