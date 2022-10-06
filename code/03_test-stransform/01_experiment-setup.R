#!/usr/bin/env R

#
# Setting up the S-transformation experiment. We wish to compare the pi_est in
# 3 conditions:
# 1. no S-transformation
# 2. "static" S-transformation (e.g. using cell type means)
# 3. "randomized" S-transformation (e.g. sampling rnorm by cell type)
# 
#

library(SingleCellExperiment)

#------
# paths
#------
proj.dpath <- "deconvo_method-paper"
source.dpath <- file.path(proj.dpath, "source")
save.dpath <- file.path(proj.dpath, "outputs/03_test-stransform")

lz.fname <- "lz-mr_expt-stransform_dlpfc-ro1.rda"
script.fnamev <- c("z_methods.R", "z_transform.R", "make_example_data.R")

sce.fpath <- "DLPFC_snRNAseq/processed-data/sce/sce_DLPFC.Rdata"

#-----
# load
#-----
# source key functions
for(scripti in script.fnamev){source(file.path(source.dpath, scripti))}

# single cell experiment
sce <- get(load(sce.fpath))

#------------------------------------
# show how it works with example data
#------------------------------------
# note:
# ltransform arg contains params for s_rescale()

# static s transform
ltrans <- list(s_rescale = list(factorv = seq(4)))
ldecon1 <- ldecon_example(k.value = 4, ltransform = ltrans)

z <- ldecon1$z_rescale
head(z)
#           k_1       k_2       k_3       k_4
#[1,] 8.7888629  8.278891 18.999112 13.543234
#[2,] 4.0212999 21.046641 15.107364 25.959136
#[3,] 8.9893978  9.059414 14.916480  2.696406
#[4,] 8.8172880  1.661184 12.278602  6.648772
#[5,] 6.2439243  1.161376 19.431195  7.825410
#[6,] 0.3801499  9.582886  9.575523 26.808364

# randomized s transform
ltrans <- list(s_rescale = list(meanv = seq(4), sdv = rep(2,4)))
ldecon2 <- ldecon_example(k.value = 4, ltransform = ltrans)
head(ldecon2$z_rescale)
#          k_1        k_2       k_3       k_4
# [1,]  0.000000 12.1385130  0.000000  0.000000
# [2,]  5.507968 41.1391628 28.127941 14.267247
# [3,] 37.536944  0.4782795  8.768361  2.898522
# [4,]  0.000000  1.7252446 20.799447 10.393403
# [5,]  5.241753  0.1964092 34.663691 12.713812
# [6,]  0.480829 14.9222000  0.000000 11.180610

#-----------------------
# conduct new experiment
#-----------------------
# params
celltype.varname <- "cellType_broad_hc"
celltype.treg.varname <- "celltype.treg"

# make new cell type -- make types per tregs paper
varv <- as.character(sce[[celltype.varname]])
varv[!varv %in% c("Excit", "Inhib", "Oligo")] <- "other"
sce[[celltype.treg.varname]] <- varv

# summarize cell type data
table(varv)
#  Astro EndoMural     Micro     Oligo       OPC     Excit     Inhib
#   3979      2157      1601     32051      1940     24809     11067
table(sce[[celltype.treg.varname]])
# Excit Inhib Oligo other
# 24809 11067 32051  9677

# get new lz for k=4
lz <- get_z_experiment(sce, ngenes.byk = 25, summary.varname = "BrNum", 
                       return.all = TRUE, type.varname = celltype.treg.varname,
                       save.dpath = save.dpath)
# inspect z.final
head(lz[["z.final"]])

# get z transformations
# z static
lz[["zs1"]] <- s_rescale(lz[["z.final"]], factorv = c(2,2,6,8))
# z randomized
lz[["zs2"]] <- s_rescale(lz[["z.final"]], meanv = c(), sdv = c())



