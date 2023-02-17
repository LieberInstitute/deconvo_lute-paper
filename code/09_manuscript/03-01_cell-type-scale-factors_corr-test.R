#!/usr/bin/env R

# Author: Sean Maden
#
# Correlation tests between cell type scale factors.
#
#

libv <- c("ggplot2", "gridExtra", "ComplexHeatmap")
sapply(libv, library, character.only = TRUE)

# manage params
handle.str <- "ro1-dlpfc"

# manage paths
code.dname <- "09_manuscript"
proj.dname <- "deconvo_method-paper"
save.dpath <- file.path(proj.dname, "outputs", code.dname)

#-------------------------
# get marker scale factors
#-------------------------

# get scale factors from snrnaseq data
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
lscale.fname <- paste0("lscale_total-counts-exprgene_",handle.str,".rda")
save(lscale, file = file.path(save.dpath, lscale.fname))

# get scale factors from halo data
halo.fname <- "halo_all.Rdata"
halo.fpath <- file.path("Human_DLPFC_Deconvolution", 
                   "processed-data", "03_HALO", halo.fname)
dfh <- get(load(halo.fpath))

#--------------------------------
# plot scale factor distributions
#--------------------------------
set.seed(0)

lscale.fname <- paste0("lscale_total-counts-exprgene_",handle.str,".rda")
lscale <- get(load(file.path(save.dpath, lscale.fname)))

# total mrna counts
variable.name <- "total.counts"
variable.str <- "total.count"
title.str <- "Total mRNA counts"
dfp <- do.call(rbind, lapply(seq(3), function(ii){
  lscale[[ii]][[variable.name]]}))
dfp$label <- paste0(dfp$type, ";", dfp$marker.type)
dfp$value <- as.numeric(dfp[,variable.str])
# downsample by label
dft <- as.data.frame(table(dfp$label))
min.cells <- min(dft[,2])
unique.labels <- unique(dfp$label)
dfp <- do.call(rbind, lapply(unique.labels, function(labeli){
  dfi <- dfp[dfp$label == labeli,]
  dfi[sample(seq(nrow(dfi)), min.cells),]
}))
# order labels
medianv <- unlist(lapply(unique.labels, function(ui){
  median(dfp[dfp$label==ui,]$value)}))
labv.order <- rev(order(medianv))
dfp$label <- factor(dfp$label, levels = unique.labels[labv.order])
# get plot object
ggvp.sn.rna <- ggplot(dfp, aes(x = label, y = value)) + theme_bw() +
  geom_violin(draw_quantiles = 0.5) + ggtitle(title.str) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggvp.sn.rna

# total expressed genes
variable.name <- "total.expr.genes"
variable.str <- "total.expr.genes"
title.str <- "Total expr. genes"
dfp <- do.call(rbind, lapply(seq(3), function(ii){
  lscale[[ii]][[variable.name]]}))
dfp$label <- paste0(dfp$type, ";", dfp$marker.type)
dfp$value <- as.numeric(dfp[,variable.str])
# downsample by label
dft <- as.data.frame(table(dfp$label))
min.cells <- min(dft[,2])
unique.labels <- unique(dfp$label)
dfp <- do.call(rbind, lapply(unique.labels, function(labeli){
  dfi <- dfp[dfp$label == labeli,]
  dfi[sample(seq(nrow(dfi)), min.cells),]
}))
# order labels
medianv <- unlist(lapply(unique.labels, function(ui){
  median(dfp[dfp$label==ui,]$value)}))
labv.order <- rev(order(medianv))
dfp$label <- factor(dfp$label, levels = unique.labels[labv.order])
# get plot object
ggvp.sn.expr <- ggplot(dfp, aes(x = label, y = value)) + theme_bw() +
  geom_violin(draw_quantiles = 0.5) + ggtitle(title.str) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggvp.sn.expr

# nucleus areas
variable.name <- "Nucleus_Area"
title.str <- paste0("Nucleus area")
dfp <- data.frame(value = dfh[,variable.name],
                  cell.type = dfh$cell_type)
dfp$label <- dfp$cell.type
# downsample by label
dft <- as.data.frame(table(dfp$label))
min.cells <- min(dft[,2])
unique.labels <- unique(dfp$label)
dfp <- do.call(rbind, lapply(unique.labels, function(labeli){
  dfi <- dfp[dfp$label == labeli,]
  dfi[sample(seq(nrow(dfi)), min.cells),]
}))
# order labels
medianv <- unlist(lapply(unique.labels, function(ui){
  median(dfp[dfp$label==ui,]$value)}))
labv.order <- rev(order(medianv))
dfp$label <- factor(dfp$label, levels = unique.labels[labv.order])
# get plot object
ggvp.img.nuc <- ggplot(dfp, aes(x = label, y = value)) + theme_bw() +
  geom_violin(draw_quantiles = 0.5) + ggtitle(title.str) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggvp.img.nuc

# nucleus areas
variable.name <- "AKT3_Copies"
title.str <- paste0("AKT3 copies")
dfp <- data.frame(value = dfh[,variable.name],
                  cell.type = dfh$cell_type)
dfp$label <- dfp$cell.type
# downsample by label
dft <- as.data.frame(table(dfp$label))
min.cells <- min(dft[,2])
unique.labels <- unique(dfp$label)
dfp <- do.call(rbind, lapply(unique.labels, function(labeli){
  dfi <- dfp[dfp$label == labeli,]
  dfi[sample(seq(nrow(dfi)), min.cells),]
}))
# order labels
medianv <- unlist(lapply(unique.labels, function(ui){
  median(dfp[dfp$label==ui,]$value)}))
labv.order <- rev(order(medianv))
dfp$label <- factor(dfp$label, levels = unique.labels[labv.order])
# get plot object
ggvp.img.akt3 <- ggplot(dfp, aes(x = label, y = value)) + theme_bw() +
  geom_violin(draw_quantiles = 0.5) + ggtitle(title.str) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggvp.img.akt3

# save composite plot
# get formatted plots list
lgg <- list(plot1 = ggvp.sn.rna + 
              theme(axis.title.x = element_blank(),
                    axis.title.y = element_blank()),
            plot2 = ggvp.sn.expr +
              theme(axis.title.x = element_blank(),
                    axis.title.y = element_blank()),
            plot3 = ggvp.img.nuc +
              theme(axis.title.x = element_blank(),
                    axis.title.y = element_blank()),
            plot4 = ggvp.img.akt3 +
              theme(axis.title.x = element_blank(),
                    axis.title.y = element_blank()))
# save jpg
jpg.fname <- "ggvp-comp_cell-scale-sn-img_ro1-dlpfc.jpg"
jpeg(file = file.path(save.dpath, jpg.fname), width = 5, height = 8,
     units = "in", res = 400)
grid.arrange(lgg[[1]], lgg[[2]], lgg[[3]], lgg[[4]],
             ncol = 1, bottom = "Label", left = "Value")
dev.off()


#---------------------
# harmonize scale data
#---------------------












