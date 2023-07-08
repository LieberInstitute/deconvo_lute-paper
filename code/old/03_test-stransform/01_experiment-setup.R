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

# saved objects
# new z datasets
lz.fname <- "lz_s-rescale-k4_dlpfc-ro1.rda"
# new pseudobulked data
lpb.fname <- "lpseudobulk_stransform-expt_dlpfc-ro1.rda"

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
save(lz, file = file.path(save.dpath, lz.fname))

# get pseudobulked tables
# realize scef subset of sce
# get scef subset
markerv <- rownames(lz$z.final)
scef <- sce[markerv,]
scef <- SingleCellExperiment(assays = list(counts = counts(scef)))
colData(scef) <- colData(sce)
# set the cell weights
datv <- c(1,2,3,3,2,3,5,2,4,2,1,1)
# get the pseudobulked data
lpb <- get_lpb(datv, scef = scef, scale.range = 100:200, 
               counts.summary.method = "mean")

head(lpb[[3]])
#               j_1       j_2       j_3
# RNF220     2.302198 1.7102273 2.2142857
# PKN2-AS1   1.824176 1.2954545 0.8452381
# MIR137HG   3.521978 1.9318182 1.7083333
# C1orf61    1.181319 1.2840909 0.9404762
# ATP1A2     1.857143 2.1590909 0.9583333
# AC011995.2 1.329670 0.9886364 0.4523810


save(lpb, file = file.path(save.dpath, lpb.fname))

#----------------------
# compare pi_est, pi_pb
#----------------------
lz <- get(load(file.path(save.dpath, lz.fname)))
lpb <- get(load(file.path(save.dpath, lpb.fname)))

# make new z tables/do s transformations
# note: order of vectors to be c(Inhib, Oligo, other, Excit)
meanv <- c(6,2,2,8); sdv <- c(2, 1, 1, 3)
# get the s-rescaled z tables
lz[["zs1"]] <- s_rescale(lz[["z.final"]], 
                         factorv = meanv)
# z randomized
lz[["zs2"]] <- s_rescale(lz[["z.final"]], 
                         meanv = meanv, sdv = sdv)

# get nnls solutions
znamev <- c("z.final", "zs1", "zs2")
y.data <- lpb[["y_data_pb"]]
# get as list
lpi <- lapply(znamev, function(znamei){
  z.data <- lz[[znamei]]
  pi.dati <- as.data.frame(do.call(rbind, lapply(seq(ncol(z.data)), function(i){
    nnls::nnls(y.data, z.data[,i])$x
  })))
  colnames(pi.dati) <- colnames(y.data)
  rownames(pi.dati) <- colnames(z.data)
  pi.dati$type <- znamei
  pi.dati
})
names(lpi) <- znamev
# get as table
df.pi <- do.call(rbind, lpi)

# get proportions from nnls
# get nnls solutions
znamev <- c("z.final", "zs1", "zs2")
y.data <- lpb[["y_data_pb"]]
# get as list
lpi <- lapply(znamev, function(znamei){
  z.data <- lz[[znamei]]
  pi.dati <- do.call(rbind, lapply(seq(ncol(z.data)), 
                                          function(i){
    nnls::nnls(y.data, z.data[,i])$x}))
  pi.dati <- apply(pi.dati, 2, function(ci){ci/sum(ci)})
  pi.dati <- as.data.frame(pi.dati)
  colnames(pi.dati) <- colnames(y.data)
  rownames(pi.dati) <- colnames(z.data)
  pi.dati$type <- znamei
  pi.dati
})
names(lpi) <- znamev
# get as table
df.pi <- do.call(rbind, lpi)

# cell types to assign
cell.typev <- c("Inhib", "Oligo", "other", "Excit")
# compare to pi_pb
pi.pb <- lpb[["pi_pb"]]
colnames(pi.pb) <- names(lpb[[2]])
rownames(pi.pb) <- cell.typev

# make new table:
# method1 cell_type sample_id pi_true pi_est diff_pb_minus_est
# method2 cell_type sample_id pi_true pi_est diff_pb_minus_est
require(reshape2)
znamev <- c("z.final", "zs1", "zs2")
y.data <- lpb[["y_data_pb"]]
pi.pb.matrix <- pi.pb
df.tall <- do.call(rbind, lapply(znamev, function(znamei){
  message(znamei)
  z.data <- lz[[znamei]]
  pi.dati <- do.call(rbind, lapply(seq(ncol(z.data)), 
                                   function(i){
                                     nnls::nnls(y.data, z.data[,i])$x}))
  pi.dati <- apply(pi.dati, 2, function(ci){ci/sum(ci)})
  colnames(pi.dati) <- colnames(y.data)
  rownames(pi.dati) <- colnames(z.data)
  # get results table
  # assign cell_type
  pi.pb.df <- as.data.frame(pi.pb.matrix)
  pi.dati <- as.data.frame(pi.dati)
  pi.dati$cell_type <- pi.pb.df$cell_type <- cell.typev
  # get tall tables
  pi.dati.tall <- melt(pi.dati, id = "cell_type")
  pi.pb.tall <- melt(pi.pb.df, id = "cell_type")
  # return tall table
  df.tall <- data.frame(cell_type = pi.dati.tall$cell_type,
                        sample_id = pi.dati.tall$variable,
                        pi_est = pi.dati.tall$value,
                        pi_true = pi.pb.tall$value,
                        pi_diff = pi.pb.tall$value-pi.dati.tall$value)
  df.tall$method <- znamei
  return(df.tall)
}))

df.tall.fname <- "df-results_s-transform-expt_dlpfc-ro1.rda"
save(df.tall, file = file.path(save.dpath, df.tall.fname))