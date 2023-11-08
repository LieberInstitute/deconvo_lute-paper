libv <- c("dplyr", "lute", "ComplexHeatmap")
sapply(libv, library, character.only=T)
csvPath<-"./data/monaco_et_al_2019/manuscript/abisseq_rnaseq_cell_types_references_k17.csv"
zTpm<-read.csv(csvPath)
zTpm <- zTpm[!duplicated(zTpm[,1]),]
rownames(zTpm) <- paste0("symbol:",zTpm[,1])
zTpm <- zTpm[,seq(2,ncol(zTpm))]
for(c in seq(ncol(zTpm))){zTpm[,c] <- as.numeric(zTpm[,c])}
zTpm <- as.matrix(zTpm)
cellSizes<-apply(as.matrix(zTpm),2,sum)
zTpmTransformed <- lute:::.zstransform(zTpm, cellSizes)
Heatmap(zTpm)
Heatmap(zTpmTransformed)
save.image(file="./env/04_transform/01_transform_script.RData")