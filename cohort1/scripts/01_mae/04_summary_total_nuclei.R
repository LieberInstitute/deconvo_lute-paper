#!/usr/bin/env R

# Author: Sean Maden
#
# Summaries of nuclei counts and proportions.
#
#

libv <- c("SingleCellExperiment", "MultiAssayExperiment")
sapply(libv, library, character.only = TRUE)

# load
load("outputs/01_mae/mae_analysis.rda")

# snrnaseq
# k4 summaries
sce <- mae[["snrnaseq.k4.all"]]
cd <- colData(sce)
dfs <- table(cd$k4, cd$Sample)
dfs <- dfs[!rownames(dfs)=="Ambiguous",] # filter
dfp <- apply(dfs, 2, function(ci){ci/sum(ci)})
apply(dfp, 1, function(ri){min(ri)})
apply(dfp, 1, function(ri){max(ri)})
apply(dfp, 1, function(ri){mean(ri)})
apply(dfp, 1, function(ri){median(ri)})

# total nuclei
#unique.samples <- unique(sc$Sample)
#mean(unlist(lapply(unique.samples, function(sample.id){
#  ncol(sc[,sc$Sample==sample.id])
#})))

#median(unlist(lapply(unique.samples, function(sample.id){
#  ncol(sc[,sc$Sample==sample.id])
#})))

#sd(unlist(lapply(unique.samples, function(sample.id){
#  ncol(sc[,sc$Sample==sample.id])
#})))