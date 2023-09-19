#!/usr/bin/env R

# Author: Sean Maden
# 
# Summary statisics for snRNAseq.
#

mae.final <- load("./outputs/01_mae/mae_allsamples.rda")
mae.sn <- list(mae.final[[1]], mae.final[[2]], mae.final[[3]])

#------------
# cell counts
#------------

lapply(c(2,3,4), function(index){
  message("K", index)
  cd <- colData(mae.sn[[index-1]])
  dft <- as.data.frame(table(cd$Sample, cd[,paste0("k", index)]))
  cell.types <- unique(cd[,paste0("k", index)])
  for(type in cell.types){
    message(type, ", mean = ", round(mean(dft[dft[,2]==type,3]), 2), "\n")
    message(type, ", median = ", round(median(dft[dft[,2]==type,3]), 2), "\n")
    message(type, ", sd = ", round(sd(dft[dft[,2]==type,3]), 2), "\n")
  }
})

#-----------------
# cell proportions
#-----------------

lapply(c(2,3,4), function(index){
  message("K", index)
  cd <- colData(mae.sn[[index-1]])
  unique.samples <- unique(cd$Sample)
  cell.types <- unique(cd[,paste0("k", index)])
  dft <- do.call(rbind, lapply(unique.samples, function(sample.id){
    cdf <- cd[cd$Sample==sample.id,]
    prop.sample.id <- prop.table(table(cdf[,paste0("k", index)]))
    prop.sample.id
  }))
  for(type in cell.types){
    message(type, ", mean = ", round(mean(dft[,type]), 2), "\n")
    message(type, ", median = ", round(median(dft[,type]), 2), "\n")
    message(type, ", sd = ", round(sd(dft[,type]), 2), "\n")
  }
})

#---------------------
# marker library sizes
#---------------------

lapply(c(2,3,4), function(index){
  message("K", index)
  sce.filt <- mae.sn[[index-1]]
  cd <- colData(sce.filt)
  unique.samples <- unique(cd$Sample)
  cell.types <- unique(cd[,paste0("k", index)])
  dft <- do.call(rbind, lapply(unique.samples, function(sample.id){
    scef <- sce.filt[,sce.filt$Sample==sample.id]
    dftf <- do.call(cbind, lapply(cell.types, function(type){
      sceff <- scef[,scef[[paste0("k",index)]]==type]
      mean(colSums(assays(sceff)[["counts"]]))
    }))
    dftf
  }))
  colnames(dft) <- cell.types
  for(type in cell.types){
    message(type, ", mean = ", round(mean(dft[,type]), 2), "\n")
    message(type, ", median = ", round(median(dft[,type]), 2), "\n")
    message(type, ", sd = ", round(sd(dft[,type]), 2), "\n")
  }
})




