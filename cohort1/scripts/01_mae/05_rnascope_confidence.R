#!/usr/bin/env R

# Author: Sean Maden
#
# Gets summaries and confidences by sample.id.
#
#
#

libv <- c("ggplot2", "reshape2")
sapply(libv, library, character.only = TRUE)

load("./outputs/01_mae/mae_allsamples_append.rda")

combo.id.variable <- "SAMPLE_ID"
confidence.variable <- "Confidence"

cd.id <- as.data.frame(colData(mae))
cd.id$confidence.circle <- cd.id$confidence.star <- "NA"

sce.img <- mae[["sce.img"]]
cd.rnascope <- colData(sce.img)
cd.rnascope$combo <- gsub(".*_", "", cd.rnascope$SAMPLE_ID)
table(cd.rnascope$combo)

cell.sizes <- metadata(mae[["snrnaseq.k2.all"]])[["cell.sizes"]]
cell.sizes <- cell.sizes[cell.sizes$ktype=="k2",]

cell.proportions <- metadata(mae[["snrnaseq.k2.all"]])[["list.df.true.k2"]]
cd.id$size.sn.glial <- cd.id$size.sn.neuron <- "NA"
cd.id$proportion.sn.glial <- cd.id$proportion.sn.neuron <- "NA"

for(sample.id in cd.id$sample.id){
  message(sample.id)
  filter.sample.rnascope <- cd.rnascope$Sample==sample.id
  cd.rnascope.iter <- cd.rnascope[filter.sample.rnascope,]
  cd.id[cd.id$sample.id==sample.id,]$confidence.circle <- cd.rnascope.iter[cd.rnascope.iter$combo=="CIRCLE",]$Confidence[1]
  cd.id[cd.id$sample.id==sample.id,]$confidence.star <- cd.rnascope.iter[cd.rnascope.iter$combo=="STAR",]$Confidence[1]
  
  cell.sizes.iter <- cell.sizes[cell.sizes$sample.id==sample.id,]
  filter.glial <- cell.sizes.iter[,2]=="glial"
  filter.neuron <- cell.sizes.iter[,2]=="neuron"
  cd.id[cd.id$sample.id==sample.id,]$size.sn.glial <- cell.sizes.iter[filter.glial,]$size[1]
  cd.id[cd.id$sample.id==sample.id,]$size.sn.neuron <- cell.sizes.iter[filter.neuron,]$size[1]
  
  if(sample.id %in% names(cell.proportions)){
    cd.id[cd.id$sample.id==sample.id,]$proportion.sn.glial <- cell.proportions[[sample.id]][["glial"]]
    cd.id[cd.id$sample.id==sample.id,]$proportion.sn.neuron <- cell.proportions[[sample.id]][["neuron"]]
  }
  
}

for(c in c(4:7)){cd.id[,c] <- as.numeric(as.character(cd.id[,c]))}

cd.id$remove.low <- cd.id$confidence.circle=="Low" | cd.id$confidence.star=="Low"

#cd.id.tall <- rbind(data.frame())
#cd.id.tall <- as.data.frame(cd.id.tall)

save(cd.id, file = "./outputs/01_mae/sample_qc_df.rda")
