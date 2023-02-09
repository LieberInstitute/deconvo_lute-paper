#!/usr/bin/env R

# Author: Sean Maden
#
# Test scale filters on abundances by cell types.
#

libv <- c("SummarizedExperiment")
sapply(libv, library, character.only = T)

# load data
lscef.fname <- "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda"
lscef <- get(load(lscef.fname))
sce <- lscef[[1]]

# get cell type abundaces by sample
cd <- colData(sce)
samp.varname <- "Sample"
celltype.varname <- "k2"
dft <- as.data.frame(table(cd[,samp.varname], cd[,celltype.varname]))

# test scale filter on lowest-abundance type (glial)
scale.thresh <- -0.9
dftf <- dft[dft$Var2=="glial",]
scalev <- scale(dftf$Freq)[,1]
which.filt <- which(scalev <= scale.thresh)
dftf[which.filt,]
# Var1  Var2 Freq
# 2 Br2720_post glial  531
# 9  Br6471_ant glial   38

# test scale filter on lowest-abundance type (glial)
scale.thresh <- -1
dftf <- dft[dft$Var2=="glial",]
scalev <- scale(dftf$Freq)[,1]
which.filt <- which(scalev <= scale.thresh)
dftf[which.filt,]
# Var1  Var2 Freq
# 9 Br6471_ant glial   38

# test scale filter on lowest-abundance type (glial)
scale.thresh <- -0.8
dftf <- dft[dft$Var2=="glial",]
scalev <- scale(dftf$Freq)[,1]
which.filt <- which(scalev <= scale.thresh)
dftf[which.filt,]
#         Var1  Var2 Freq
# 2 Br2720_post glial  531
# 4  Br2743_mid glial  572
# 9  Br6471_ant glial   38
