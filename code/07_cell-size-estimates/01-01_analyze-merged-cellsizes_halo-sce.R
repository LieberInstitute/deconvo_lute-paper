#!/usr/bin/env R

# Author: Sean Maden
#

#----------
# load data
#----------
# read output df objects
read.fname <- "df-merge-cellsize_halo-sce.rda"
read.dpath <- out.dpath <- file.path("deconvo_method-paper", "outputs", 
                                     "07_cell-size-estimates")
dfm <- get(load(file.path(read.dpath, read.fname)))

#----------------------
# analysis by cell type
#----------------------
ctv <- c("Excit", "Inhib", "Oligo")
cti <- "Excit"

lcor <- lapply(ctv, function(cti){
  dffm <- dfm[,grepl(cti, colnames(dfm))]
  for(c in seq(ncol(dffm))){dffm[,c] <- as.numeric(dffm[,c])}
  cor(dffm, use = "pairwise.complete.obs", method = "spearman")
})


