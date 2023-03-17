#!/usr/bin/env R

# Author: Sean Maden
#
# Pre-filter and process bulk RNA-seq data.
#

libv <- c("SummarizedExperiment")
sapply(libv, library, character.only = T)

#----------
# load data
#----------
# load bulk data
rse.filename <- "rse_gene.Rdata"
load.path <- file.path("Human_DLPFC_Deconvolution",
                       "processed-data",
                       "01_SPEAQeasy",
                       "round2_v40_2022-07-06",
                       "rse")
rse <- get(load(file.path(path, rse.filename)))

# get save directory path
save.path <- file.path("deconvo_method-paper", "outputs", 
                       "11_bulk-expression-analysis")

#-------------
# filter genes
#-------------
# set type labels to filter
types.protein <- c("protein_coding")
types.nonpolya <- c("lncRNA", "Mt_rRNA", "rRNA", "Mt_tRNA")
types.include <- c(types.protein, types.nonpolya)
# get gene types
rd <- rowData(rse)
types.vector <- rd$gene_type

#-----------------------
# get rescale transforms
#-----------------------

#-----------------------
# save preprocessed data
#-----------------------
rse.filename <- "rse_processed.rda"
save(rse, file = file.path(save.path, rse.filename))

# filter protein-coding loci
## filter protein-coding
## want to include genes that aren't polyA selected, e.g. MT genes, ribosomal RNA, histones, long NC RNAs ...
## check for MT gene significance in global DEGs

# consider other uses for non-coding loci, discuss with KM
## filter miRNA
## filter longncRNA
## filter mtRNA

# note scale from tx decon paper -- TPM? Also note other methods details we can apply
## Does SPEQeasy use STAR under the hood, similar to this paper?

# add normalized counts?
# add logcounts?

# notes for deg tests
## are markers among within-bulk degs?

# review paper Hippen et al 2023 -- tx decon by bulk processing protocol