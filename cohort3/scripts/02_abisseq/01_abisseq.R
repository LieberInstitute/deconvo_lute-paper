#!/usr/bin/env R

# Author: Sean Maden
#
# Performs deconvolution with the ABIS-seq signature matrix reference.
#
#
#

libv <- c("lute")
sapply(libv, library, character.only = TRUE)

#-----
# load
#-----
zref <- read.csv("./data/zref_abisseq.csv")
ptrue <- read.csv("./data/ptrue_abis.csv")
load("./env/01_tpm_summaries/01_read_script.RData")

#--------------
# format inputs
#--------------
# zref
rownames(zref) <- zref[,1]
zref <- zref[,c(2:ncol(zref))]

#------------
# run deconvo
#------------

list.result <- lute(
  z = as.matrix(zref), 
  y = assays(se)[["tpm"]], 
  typemarker.algorithm = NULL
  )




#-----
# save
#-----
save.image("./env/02_abisseq/01_abisseq_script.RData")
