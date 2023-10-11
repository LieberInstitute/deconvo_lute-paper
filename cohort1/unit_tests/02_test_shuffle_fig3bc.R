#!/usr/bin/env R
# Author: Sean Maden

expected.rnascope.sample.count <- 1
expected.snrnaseq.sample.count <- 1

load("./env/03_shuffle/00_fig3bc_script.RData")

#-----------------------------------
# test 1, number of RNAscope samples
#-----------------------------------
# dim(dfs.rn)
# length(unique(dfs.rn$Sample))

length(unique(dfs.rn$sample.id))==expected.rnascope.sample.count

#-----------------------------------
# test 2, number of RNAscope samples
#-----------------------------------
# dim(sce.k2)

length(sce.k2$Sample)==expected.snrnaseq.sample.count
