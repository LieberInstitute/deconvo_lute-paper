#!/usr/bin/env R

# Author: Sean Maden
#
# Gets heatmap data before and after cell type scale factor transformations.
#

libv <- c("dplyr", "lute", "ComplexHeatmap")
sapply(libv, library, character.only=T)
# param
markersPerType<-20
# load
csvPath<-"./data/monaco_et_al_2019/manuscript/abisseq_rnaseq_cell_types_references_k17.csv"
refTpm<-read.csv(csvPath)
# format
refTpm <- refTpm[!duplicated(refTpm[,1]),]
rownames(refTpm)<-refTpm[,1]
refTpm <- refTpm[,seq(2,ncol(refTpm))]
for(c in seq(ncol(refTpm))){refTpm[,c] <- as.numeric(refTpm[,c])}

# TPM heatmap data
refTpm <- as.matrix(refTpm)
# transform
cellSizes<-apply(as.matrix(refTpm),2,sum)
refTpmTransformed <- lute:::.zstransform(refTpm, cellSizes)

# log2 TPM + 1 heatmap data
refTpmLog2 <- as.data.frame(t(apply(refTpm,1,function(ci){
  log2(ci+1)
})))
colnames(refTpmLog2) <- colnames(refTpm)
rownames(refTpmLog2) <- rownames(refTpm)
# transform
cellSizesLog2Tpm<-apply(as.matrix(refTpm),2,sum)
refTpmLog2Transformed <- lute:::.zstransform(refTpmLog2, cellSizesLog2Tpm)

# filter markers
# TEST: does listMarkers change with TPM vs log2TPM+1?
listMarkers <-lapply(
  colnames(refTpm), function(cellType){
    cellTypeTrue<-colnames(refTpm)==cellType
    cellTypeFalse<-!colnames(refTpm)==cellType
    tpmTrue<-refTpm[,cellTypeTrue]
    tpmFalse<-apply(refTpm[,cellTypeFalse],1,median)
    tpmDiff<-tpmTrue-tpmFalse
    tpmAbsDiff<-abs(tpmDiff)
    tpmAbsDiffOrdered<-tpmAbsDiff[rev(order(tpmAbsDiff))]
    return(names(tpmAbsDiffOrdered[seq(markersPerType)]))
  }
)
names(listMarkers)<-colnames(refTpm)

# tables for heatmaps
refTpmFilter<-refTpm[unlist(listMarkers),]
refTpmTransformedFilter<-refTpmTransformed[unlist(listMarkers),]

# heatmaps
Heatmap(refTpm)
Heatmap(scale(refTpm))

# save
save.image(file="./env/04_transform/01_transform_script.RData")