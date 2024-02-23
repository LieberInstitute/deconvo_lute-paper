#!/usr/bin/env R

# Author: Sean Maden
#
# Gets heatmap data for top 10 markers per cell type in ABIS-seq reference.
#

libv <- c("dplyr", "lute", "ComplexHeatmap")
sapply(libv, library, character.only=T)

# param
markersPerType <- 10

# load
load("./env/01_pseudobulk/01_read_script.RData")

# format
#refTpm <- refTpm[!duplicated(refTpm[,1]),]
#for(c in seq(ncol(refTpm))){refTpm[,c] <- as.numeric(refTpm[,c])}

# TPM heatmap data
refTpm <- experimentData$z
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
cellSizesLog2Tpm<-apply(as.matrix(refTpmLog2),2,sum)
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
refTpmLog2Filter <- refTpmLog2[unlist(listMarkers),]
refTpmLog2TransformedFilter <- refTpmLog2Transformed[unlist(listMarkers),]

# heatmaps
Heatmap(refTpm)
Heatmap(scale(refTpm))

# save
save.image(file="./env/01_pseudobulk/02_markers_script.RData")