#!/usr/bin/env R

# 
#
#

libv <- c("ggplot2", "reshape2")
sapply(libv, library, character.only = TRUE)

load("./outputs/01_mae/mae_allsamples_append.rda")

combo.id.variable <- "SAMPLE_ID"
confidence.variable <- "Confidence"

cd.id <- colData(mae)
sce.img <- mae[["sce.img"]]
cd.rnascope <- colData(sce.img)
cd$confidence.circle <- cd$confidence.star <- "NA"

sce.img$combo <- gsub(".*_", "", sce.img[[combo.id.variable]]) 

for(sample.id in cd$sample.id){
  cd.id[cd.id]$confidence <- colData(sce.img)[sce.img$combo=="STAR",confidence.variable][1]
  
}



