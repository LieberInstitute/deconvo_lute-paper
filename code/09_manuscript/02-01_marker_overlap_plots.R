#!/usr/bin/env R

#
# Marker overlap heatmap annotation and upset plots.
#
#

libv <- c("scuttle", "ComplexHeatmap", "UpSetR")
sapply(libv, library, character.only = TRUE)

set.seed(0) # seed for random colors

#----------
# load data
#----------
# ro1 dlpfc sce, markers
fname1 <- "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda"
lscef1 <- get(load(fname1))

# mrb dlpfc sce, markers
fname2 <- "list-scef_markers-k2-k3-k4_mrb-dlpfc.rda"
lscef2 <- get(load(fname2))

# lit annotations
llit.fnamev <- c("lctmarkers-gsheet-dlpfc.rda")
llit <- get(load(file.path(llit.fnamev)))

#-------------------
# format data labels
#-------------------
# format mrb labels
# rename "other" -> "glial"/"non_oligo_glial"
markerv <- names(lscef2)
for(mi in markerv){
  cd <- colData(lscef2[[mi]])
  filt <- cd[,mi] == "other"
  if(mi == "k3"){
    mstr <- "non_oligo_glial"
  } else{
    mstr <- "glial"
  }
  cd[filt, mi] <- mstr
  colData(lscef2[[mi]]) <- cd
}

#-------------
# overlap data
#-------------

markerv <- names(lscef1)

# get full gene set
all.genes <- unique(unlist(lapply(mv, function(mi){
  rownames(lscef1[[mi]])})))

# get full dfol
dfol1 <- data.frame(gene.name = all.genes)
dfol2 <- do.call(cbind, lapply(mv, function(mi){
  dfi <- dfol1; dfi[,mi] <- "NA"
  sce <- lscef1[[mi]]
  sce <- logNormCounts(sce)
  sce.means <- aggregateAcrossCells(sce, ids = sce[[mi]],
                                    use.assay.type = "logcounts")
  me <- assays(sce.means)$logcounts
  cnv <- colnames(me)
  marker.labels <- unlist(apply(me, 1, function(ri){
    cnv[which(ri==max(ri))]}))
  for(genei in names(marker.labels)){
    filt <- which(dfi[,1]==genei)
    dfi[filt, mi] <- marker.labels[genei]
  }
  return(dfi[,2])
}))
dfol3 <- do.call(cbind, lapply(mv, function(mi){
  dfi <- dfol1; 
  cname <- paste0(mi, "_mrb")
  dfi[, cname] <- "NA"
  sce <- lscef2[[mi]]
  sce <- logNormCounts(sce)
  sce.means <- aggregateAcrossCells(sce, ids = sce[[mi]],
                                    use.assay.type = "logcounts")
  me <- assays(sce.means)$logcounts
  cnv <- colnames(me)
  marker.labels <- unlist(apply(me, 1, function(ri){
    cnv[which(ri==max(ri))]}))
  for(genei in names(marker.labels)){
    filt <- which(dfi[,1]==genei)
    dfi[filt, cname] <- marker.labels[genei]
  }
  return(dfi[,2])
}))
dfol <- cbind(dfol1, dfol2, dfol3)
colnames(dfol) <- c("gene.name", "k2", "k3", "k4", 
                    "k2_mrb", "k3_mrb", "k4_mrb")

dfol[,2] <- factor(dfol[,2], levels = c("neuron", "glial", "NA"))
dfol[,3] <- factor(dfol[,3], levels = c("Excit", "Inhib", "glial", "NA"))
dfol[,4] <- factor(dfol[,4], levels = c("Excit", "Inhib", "Oligo", 
                                        "non_oligo_glial", "NA"))

dfol[,5] <- factor(dfol[,5], levels = c("neuron", "glial", "NA"))
dfol[,6] <- factor(dfol[,6], levels = c("Excit", "Inhib", "glial", "NA"))
dfol[,7] <- factor(dfol[,7], levels = c("Excit", "Inhib", "Oligo", 
                                        "non_oligo_glial", "NA"))

order.dfol <- order(dfol[,2], 
                    dfol[,3], 
                    dfol[,4], 
                    dfol[,5])
dfol <- dfol[order.dfol,]

cnv <- colnames(dfol)[2:ncol(dfol)]
catv <- unique(unlist(dfol[,cnv]))
colv <- c("blue", "orange", "gray", 
          "green", "yellow", "purple", 
          "brown")
colstr <- paste0(catv, " = ", colv)
collist <- lapply(cnv, function(ci){
  catvf <- unique(dfol[,ci])
  which.catvf <- which(catv %in% catvf)
  colvf <- colv[which.catvf]
  catvf <- catv[which.catvf]
  eval.str <- paste0("c(", 
                     paste0("'",catvf, "'='", 
                            colvf, "'", collapse = ","), 
                     ")")
  eval(parse(text = eval.str))
})
names(collist) <- cnv

# make horizontal annotation
ann = HeatmapAnnotation(k2 = dfol[,2],
                        k3 = dfol[,3],
                        k4 = dfol[,4],
                        k2_mrb = dfol[,5],
                        k3_mrb = dfol[,6],
                        k4_mrb = dfol[,7],
                        foo = anno_text(dfol[,1], 
                                        gp = gpar(fontsize = 10)),
                        col = collist,
                        gp = gpar(col = "black"),
                        annotation_name_side = "left")
plot(ann, heatmap_legend_side = "top")

# make vertical annotation
ann = rowAnnotation(foo = anno_text(dfol[,1], 
                                    gp = gpar(fontsize = 5),
                                    just = "center"),
                    k2 = dfol[,2],
                        k3 = dfol[,3],
                        k4 = dfol[,4],
                        k2_mrb = dfol[,5],
                        k3_mrb = dfol[,6],
                        k4_mrb = dfol[,7],
                        
                        col = collist,
                        gp = gpar(col = "black"),
                    annotation_name_side = "top",
                    annotation_name_rot = 45)
plot(ann, heatmap_legend_side = "right")

# only annotation for k2
dfol.filt <- dfol[!dfol[,2]=="NA",]
# make vertical annotation
ann.filt = rowAnnotation(foo = anno_text(dfol.filt[,1], 
                                    gp = gpar(fontsize = 12),
                                    just = "center",
                                    location = 0.5),
                    k2 = dfol.filt[,2],
                    k3 = dfol.filt[,3],
                    k4 = dfol.filt[,4],
                    k2_mrb = dfol.filt[,5],
                    k3_mrb = dfol.filt[,6],
                    k4_mrb = dfol.filt[,7],
                    col = collist,
                    gp = gpar(col = "black"),
                    annotation_name_side = "top",
                    annotation_name_rot = 45)
plot(ann.filt, heatmap_legend_side = "right")

# save new plot
fname <- "hmblock-overlaps_k2-mrb-ro1.jpg"
jpeg(fname, file = file.path(fname), 
     width = 6, height = 7.5,
     units = "in", res = 400)
plot(ann.filt, heatmap_legend_side = "right")
dev.off()


#------------
# upset plots
#------------
# upset overall
lgene1 <- lapply(lscef1, function(ii){rownames(ii)})
lgene2 <- lapply(lscef2, function(ii){rownames(ii)})
names(lgene1) <- paste0(names(lscef1))
names(lgene2) <- paste0(names(lscef2), "_mrb")
lupset <- do.call(c, list(lgene1, lgene2))
# save new upset plot
fname <- "upset-all-markers_k2-3-4_mrb-ro1-dlpfc.jpg"
jpeg(file = fname, width = 6.5, height = 2.8, units = "in", res = 400)
upset(fromList(lupset), order.by = "freq", nsets = length(lupset),
      number.angles = 45)
dev.off()

#----------------------------
# append lit overlaps to dfol
#----------------------------

ldat <- llit$`Oord and Aberg 2022`
genev <- unlist(ldat)
oord <- unlist(lapply(dfol$gene.name, function(genei){
  names(genev[genei])
}))

ldat <- llit$`Sjostedt et al 2015`
genev <- unlist(ldat)
sjo <- unlist(lapply(dfol$gene.name, function(genei){
  unique(names(genev[genei]))
}))

ldat <- llit$`Mathys et al 2019`
genev <- unlist(ldat)
mty <- unlist(lapply(dfol$gene.name, function(genei){
  unique(names(genev[genei]))
}))

ldat <- llit$`Maynard et al 2021`
genev <- unlist(ldat)
myd <- unlist(lapply(dfol$gene.name, function(genei){
  unique(names(genev[genei]))
}))