#!/usr/bin/env R

# Author: Sean Maden
#
# Correlation tests between cell type scale factors.
#
#

# manage paths
code.dname <- "09_manuscript"
proj.dname <- "deconvo_method-paper"
save.dpath <- file.path(proj.dname, "outputs", code.dname)

#-------------------------
# get marker scale factors
#-------------------------
handle.str <- "ro1-dlpfc"
marker.typev <- c("k2", "k3", "k4")

# get scale factors
lscale <- lapply(marker.typev, function(markeri){
  message("loading the data...")
  sce.fname <- paste0("sce_marker-adj-",markeri,"_",
                      handle.str,".rda")
  sce.fpath <- file.path(save.dpath, sce.fname)
  scei <- get(load(sce.fpath))
  
  # get cell size scale factors
  type.vector <- scei[[markeri]]
  unique.types <- unique(type.vector)
  # get tall table of total mrna by type
  df.tc <- do.call(rbind, lapply(unique.types, function(typei){
    ctf <- counts(scei[,scei[[markeri]]==typei])
    dfi <- data.frame(total.count = as.character(unlist(colSums(ctf))))
    dfi$type <- typei; dfi$marker.type <- markeri; return(dfi)
  }))
  
  # get tall table of total expressed genes by type
  df.eg <- do.call(rbind, lapply(unique.types, function(typei){
    ctf <- counts(scei[,scei[[markeri]]==typei])
    m.eg <- apply(ctf,2,function(ci){length(ci[ci>0])})
    dfi <- data.frame(total.expr.genes = as.numeric(unlist(m.eg)))
    dfi$type <- typei; dfi$marker.type <- markeri; return(dfi)
  }))
  
  lr <- list(total.counts = df.tc, total.expr.genes = df.eg)
  return(lr)
})
names(lscale) <- marker.typev

# save
fname <- paste0("lscale_total-counts-exprgene_",handle.str,".rda")
save(lscale, file = file.path(save.dpath, fname))