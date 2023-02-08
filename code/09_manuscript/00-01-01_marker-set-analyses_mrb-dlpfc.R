#!/usr/bin/env R

#
#
#
#

handle.str <- "mrb-dlpfc"

libv <- c('lute', "ggplot2", "SingleCellExperiment", "limma",
          "SummarizedExperiment")
sapply(libv, library, character.only = TRUE)

#----------
# load data
#----------
# get save dpath
code.dname <- "09_manuscript"
proj.dname <- "deconvo_method-paper"
save.dpath <- file.path(proj.dname, "outputs", code.dname)

#----------------------
# get marker expression
#----------------------
marker.typev <- c("k2", "k3", "k4")

lscef <- lapply(marker.typev, function(markeri){
  message("loading the data...")
  sce.fname <- paste0("sce_marker-adj-",markeri,"_",
                      handle.str,".rda")
  sce.fpath <- file.path(save.dpath, sce.fname)
  scei <- get(load(sce.fpath))
  
  # get marker expr
  mr <- metadata(sce)[["markers"]][["top"]] # top markers
  scef <- scei[mr$gene,]
  return(scef)
})
names(lscef) <- marker.typev

# save
fname <- paste0("list-scef_markers-k2-k3-k4_",handle.str,".rda")
save(lscef, file = file.path(save.dpath, fname))

#------------------------------
# heatmaps of scaled expression
#------------------------------
library(ComplexHeatmap)

set.seed(0)

# params
batchvar <- "donor"
celltypevar <- "k2"
assayname.adj <- "counts_adj"
hm.name <- "Adj. counts,\nscaled"
hm.height <- 5
hm.width <- 6

for(ii in seq(length(lscef))){
  markeri <- names(lscef)[ii]
  message("working on marker ", markeri, "...")
  scei <- lscef[[ii]]
  # iter
  # scei <- lscef[[1]]
  # get marker types
  scef.means2 <- aggregateAcrossCells(x = scei, id = scei[[celltypevar]],
                                      use.assay.type = assayname.adj)
  me2 <- assays(scef.means2)[[assayname.adj]]; cnv <- colnames(me2)
  markerv <- unlist(apply(me2, 1, function(ri){cnv[which(ri==max(ri))]}))
  markerv <- as.character(markerv)
  
  # plot get means by donors
  agg.var <- paste0(scei[[batchvar]], ";", scei[[celltypevar]])
  scef.means <- aggregateAcrossCells(x = scei, id = agg.var, 
                                     use.assay.type = assayname.adj)
  hm <- assays(scef.means)[[assayname.adj]]
  hm <- scale(hm)
  topanno <- HeatmapAnnotation(group = gsub(";.*", "", colnames(hm)), 
                               celltype = gsub(".*;", "", colnames(hm)))
  leftanno <- rowAnnotation(marker_type = markerv)
  hm <- Heatmap(matrix = hm, name = hm.name, show_column_dend = FALSE, 
                top_annotation = topanno, left_annotation = leftanno)
  # save new plot
  fname <- paste0("hm-markers-scaled_by-celltype-donor_marker-", 
                  markeri,"_mrb-dlpfc.jpg")
  jpeg(file.path(save.dpath, fname), 
       width = hm.width, height = hm.height, 
       units = "in", res = 400)
  print(hm); dev.off()
  
  # plot get means by type
  agg.var <- scei[[celltypevar]]
  scef.means <- aggregateAcrossCells(x = scei, id = agg.var, 
                                     use.assay.type = assayname.adj)
  hm <- assays(scef.means)[[assayname.adj]]
  hm <- scale(hm)
  topanno <- HeatmapAnnotation(celltype = gsub(".*;", "", colnames(hm)))
  leftanno <- rowAnnotation(marker_type = markerv)
  hm <- Heatmap(matrix = hm, name = hm.name, show_column_dend = FALSE, 
                top_annotation = topanno, left_annotation = leftanno)
  
  # save new plot
  fname <- paste0("hm-markers-scaled_by-celltype_marker-", 
                  markeri,"_mrb-dlpfc.jpg")
  jpeg(file.path(save.dpath, fname), 
       width = hm.width, height = hm.height, 
       units = "in", res = 400)
  print(hm); dev.off()
}

#------------
# upset plots
#------------
library(UpSetR)

# upset -- overall
lmarker <- lapply(lscef, function(scei){rownames(scei)})
names(lmarker) <- names(lscef)
m <- fromList(lmarker); class(m)
m <- as.matrix(m)
upset(fromList(lmarker), order.by = "freq")

# upset -- by cell type -- neuron

#----------------------------------------------
# enrichment and overlap of independent markers
#----------------------------------------------
# load tables
llit.fnamev <- c("lctmarkers-gsheet-dlpfc.rda")
llit <- get(load(file.path(save.dpath, llit.fnamev)))

# show cell types in lit marker list
for(ii in seq(length(llit))){
  message(ii, ", ", names(llit)[ii], ": ", 
          paste0(names(llit[[ii]]), collapse = "; "))}
# 1, Faraco et al 2017: MICRO
# 2, Swanson et al 2020: MICRO
# 3, Swanson et al 2022: MICRO
# 4, De Picker and Morrens 2020: MICRO
# 5, Notter et al 2020: MICRO
# 6, Webster et al 2005: ASTRO
# 7, Murphy et al 2020: ASTRO
# 8, Tocker et al 2018: ASTRO
# 9, Sjostedt et al 2015: ASTRO; EXCIT; NEUR; OLIGO; MICRO; ENDO
# 10, Roberts et al 2004: MACRO
# 11, Lisi et al 2017: MACRO
# 12, Zhu et al 2022: MACRO
# 13, Maynard et al 2021: INHIB
# 14, Oord and Aberg 2022: OLIGO; EXCIT
# 15, Nagy et al 2020: OPC
# 16, Mathys et al 2019: MICRO; OPC; ENDO; PERI

# upset plot -- k2 neurons
lmi <- list()
lmi[["k2"]] <- lmarker$k2
lmi[["oord"]] <- llit[[14]]$EXCIT
lmi[["sjostedt"]] <- unique(c(llit[[9]]$EXCIT, llit[[9]]$NEUR))
lmi[["maynard"]] <- unique(c(llit[[13]]$INHIB))
# plot
upset(fromList(lmi), order.by = "freq")

# upset plot -- k2 other
lmi <- list()
lmi[["k2"]] <- lmarker$k2
lmi[["oord"]] <- llit[[14]]$OLIGO
lmi[["sjostedt"]] <- unique(c(llit[[9]]$ASTRO, 
                              llit[[9]]$OLIGO, 
                              llit[[9]]$MICRO, 
                              llit[[9]]$ENDO))
lmi[["mathys"]] <- unique(c(llit[[16]]$MICRO,
                            llit[[16]]$OPC,
                            llit[[16]]$ENDO,
                            llit[[16]]$PERI))
# plot
upset(fromList(lmi), order.by = "freq")

# upset plot -- k3 excit
lmi <- list()
lmi[["k3"]] <- lmarker$k3
lmi[["oord"]] <- llit[[14]]$EXCIT
lmi[["sjostedt"]] <- unique(c(llit[[9]]$EXCIT, llit[[9]]$NEUR))
# plot
upset(fromList(lmi), order.by = "freq")

# upset plot -- k3 inhib
lmi <- list()
lmi[["k3"]] <- lmarker$k3
lmi[["sjostedt"]] <- unique(c(llit[[9]]$NEUR))
lmi[["maynard"]] <- unique(c(llit[[13]]$INHIB))
# plot
upset(fromList(lmi), order.by = "freq")

# upset plot -- k3 other
lmi <- list()
lmi[["k3"]] <- lmarker$k3
lmi[["oord"]] <- llit[[14]]$OLIGO
lmi[["sjostedt"]] <- unique(c(llit[[9]]$ASTRO, 
                              llit[[9]]$OLIGO, 
                              llit[[9]]$MICRO, 
                              llit[[9]]$ENDO))
lmi[["mathys"]] <- unique(c(llit[[16]]$MICRO,
                            llit[[16]]$OPC,
                            llit[[16]]$ENDO,
                            llit[[16]]$PERI))
# plot
upset(fromList(lmi), order.by = "freq")

# upset plot -- k4 excit
lmi <- list()
lmi[["k4"]] <- lmarker$k4
lmi[["oord"]] <- llit[[14]]$EXCIT
lmi[["sjostedt"]] <- unique(c(llit[[9]]$EXCIT, llit[[9]]$NEUR))
# plot
upset(fromList(lmi), order.by = "freq")

# upset plot -- k4 excit
lmi <- list()
lmi[["k4"]] <- lmarker$k4
lmi[["sjostedt"]] <- unique(c(llit[[9]]$NEUR))
lmi[["maynard"]] <- unique(c(llit[[13]]$INHIB))
# plot
upset(fromList(lmi), order.by = "freq")

# upset plot -- k4 oligo
lmi <- list()
lmi[["k4"]] <- lmarker$k4
lmi[["oord"]] <- llit[[14]]$OLIGO
lmi[["sjostedt"]] <- unique(c(llit[[9]]$OLIGO))
# plot
upset(fromList(lmi), order.by = "freq")

# upset plot -- k4 other
lmi <- list()
lmi[["k4"]] <- lmarker$k4
lmi[["sjostedt"]] <- unique(c(llit[[9]]$ASTRO,
                              llit[[9]]$MICRO, 
                              llit[[9]]$ENDO))
lmi[["mathys"]] <- unique(c(llit[[16]]$MICRO,
                            llit[[16]]$OPC,
                            llit[[16]]$ENDO,
                            llit[[16]]$PERI))
# plot
upset(fromList(lmi), order.by = "freq")