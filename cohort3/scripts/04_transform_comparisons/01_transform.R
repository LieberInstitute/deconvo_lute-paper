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
refTpm <- as.matrix(refTpm)
# transform
cellSizes<-apply(as.matrix(refTpm),2,sum)
refTpmTransformed <- lute:::.zstransform(refTpm, cellSizes)
# filter
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
Heatmap(refTpmTransformed)
save.image(file="./env/04_transform/01_transform_script.RData")