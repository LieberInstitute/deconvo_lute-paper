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

sce.fpath <- "DLPFC_snRNAseq/processed-data/sce/sce_DLPFC.Rdata"

script.fnamev <- c("z_methods.R", "z_transform.R", "make_example_data.R", 
                   "z_figures.R", "y_methods.R")

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
dim(lz[[3]]) # [1] 100   4
head(lz[[3]])
#                Inhib      Oligo      other    Excit
# RNF220     0.6776153 6.65279170 0.88312716 1.882595
# PKN2-AS1   0.6086981 0.12758524 0.14509842 3.043970
# MIR137HG   1.3667413 0.21631744 0.18032267 5.961459
# C1orf61    1.0006925 0.32648768 2.18928067 1.026195
# ATP1A2     0.1820084 0.38670884 4.05657646 0.407281
# AC011995.2 0.2868810 0.09214082 0.03466446 2.787617

# save new lz
save(lz, file = file.path(save.dpath, "lz_s-rescale-k4_dlpfc-ro1.rda"))

# get z transformations
# note: order of vectors to be c(Inhib, Oligo, other, Excit)
meanv <- c(6,2,2,8); sdv <- c(2, 1, 1, 3)
# z static
lz[["zs1"]] <- s_rescale(lz[["z.final"]], factorv = meanv)
# z randomized
lz[["zs2"]] <- s_rescale(lz[["z.final"]], meanv = meanv, sdv = sdv)

# get pseudobulked tables
markerv <- rownames(lz$z.final)
datv <- c(1,2,3,3,2,3,5,2,4,2,1,1)

# realize scef subset of sce
scef <- sce[markerv,]
scef <- SingleCellExperiment(assays = list(counts = counts(scef)))
colData(scef) <- colData(sce)
lpb <- get_lpb(datv, scef = scef, scale.range = 100:200, 
               counts.summary.method = "mean")

DelayedMatrixStats::rowMeans2(lpb[[1]][[1]])

# test duration of process for small table
t1 <- Sys.time()
mm <- DelayedArray::rowMeans(lpb[[1]][[1]])
t2 <- Sys.time()
t1-t2

rowMeans(as.matrix(lpb[[1]][[1]]))

names(lpb) # [1] "listed_counts_pb" "y_data_pb"        "pi_pb"
lpb$pi_pb
#             j_1       j_2   j_3
# Inhib 0.1111111 0.1666667 0.500
# Oligo 0.2222222 0.2500000 0.250
# other 0.3333333 0.4166667 0.125
# Excit 0.3333333 0.1666667 0.125

# get pi_est series
y.data <- lpb[["y_data_pb"]]

pi.dat <- do.call(rbind, lapply(seq(ncol(z.data)), function(i){
  nnls::nnls(y.data, z.data[,i])$x
}))



# eval pi differences

