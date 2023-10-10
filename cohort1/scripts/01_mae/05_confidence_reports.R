cd.id <- colData(mae)
sce.img <- mae[["sce.img"]]

cd.rnascope <- colData(sce.img)
cd$confidence.circle <- cd$confidence.star <- "NA"
combo.id.variable <- "SAMPLE_ID"
confidence.variable <- "Confidence"
sce.img$combo <- gsub(".*_", "", sce.img$SAMPLE_ID) 

for(sample.id in cd$sample.id){
  cd.id[cd.id]$confidence <- colData(sce.img)[sce.img$combo=="STAR",confidence.variable][1]
  
}