#!/usr/bin/env R

#
# Make marker heatmaps for presentation
#
#

libv <- c("scuttle", "ComplexHeatmap")
sapply(libv, library, character.only = TRUE)

set.seed(0) # seed for random colors

#--------------
# load data
#-------------- 
fname1 <- "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda"
lscef1 <- get(load(fname1))

fname2 <- "list-scef_markers-k2-k3-k4_mrb-dlpfc.rda"
lscef2 <- get(load(fname2))

#------------------------------------
# heatmaps for ro1 -- sample;celltype
#------------------------------------
proj.str <- "ro1-dlpfc"

# params
group.variable <- "Sample"
celltype.variable <- "k2"
assayname.adj <- "counts_adj"
hm.name <- "Scaled\nlogcounts"
hm.height <- 5
hm.width <- 6

# get heatmap data
sce <- lscef1[[celltype.variable]]
sce <- logNormCounts(sce, assay.type = assayname.adj)
celltype.vector <- sce[[celltype.variable]]
group.vector <- paste0(sce[[group.variable]], ";", celltype.vector)
sce.means <- aggregateAcrossCells(x = sce, ids = group.vector,
                                  use.assay.type = "logcounts")
hm1 <- assays(sce.means)[["logcounts"]]
hm1 <- scale(hm1)

# get top anno
cnv <- colnames(sce.means)
sampv <- gsub(";.*", "", cnv)
brv <- gsub("_.*", "", sampv)
locv <- gsub(".*_", "", sampv)
ctv <- gsub(".*;", "", cnv)
collist <- list(brnum = c("Br2720" = "red",
                           "Br2743" = "purple",
                           "Br3942" = "green",
                           "Br6423" = "blue",
                           "Br6471" = "pink",
                           "Br8325" = "orange",
                           "Br8492" = "gray"),
                location = c("ant" = "brown",
                             "mid" = "green",
                             "post" = "orange"),
                celltype = c("glial" = "yellow",
                             "neuron" = "blue"))
topanno <- HeatmapAnnotation(brnum = brv, 
                             location = locv,
                             celltype = ctv,
                             col = collist,
                             gp = gpar(col = "black"))

# get left anno
# get marker status
sce <- logNormCounts(sce, assay.type = assayname.adj)
sce.means2 <- aggregateAcrossCells(x = sce, ids = celltype.vector,
                                   use.assay.type = "logcounts")
me2 <- assays(sce.means2)[["logcounts"]]
cnv <- colnames(me2)
marker.vector <- unlist(apply(me2, 1, function(ri){
  cnv[which(ri==max(ri))]}))
marker.vector <- as.character(marker.vector)
# get overlapping markers -- k2
mv1 <- rownames(lscef1[[1]])
mv2 <- rownames(lscef2[[1]])
overlapping.markers <- intersect(mv1, mv2)
ol.variable <- ifelse(rownames(hm1) %in% overlapping.markers, 
                      "true", "false")
leftanno <- rowAnnotation(celltype = marker.vector,
                          overlapping = ol.variable,
                          col = list(celltype = c("neuron" = "blue",
                                                  "glial" = "yellow"),
                                     overlapping = c("true" = "purple",
                                                     "false" = "orange")),
                          show_legend = TRUE,
                          gp = gpar(col = "black"))

# heatmap --- sample;celltype
hm.sample.celltype <- Heatmap(matrix = hm1, 
              name = hm.name, 
              show_column_dend = FALSE, 
              show_row_dend = FALSE,
              top_annotation = topanno, 
              left_annotation = leftanno,
              show_column_names = FALSE,
              show_row_names = TRUE,
              rect_gp = gpar(col = "black", lwd = 1))
# save new plot
fname <- paste0("hm-markers_sample-celltype_marker-", 
                celltype.variable,"_", proj.str,".jpg")
jpeg(fname, width = 6, 
     height = 8.5, units = "in", 
     res = 400)
print(hm.sample.celltype); dev.off()


# heatmap -- celltype
# get heatmap data
hm2 <- scale(me2)
# get top anno
ctv <- colnames(hm2)
collist <- list(celltype = c("glial" = "yellow", "neuron" = "blue"))
topanno <- HeatmapAnnotation(celltype = ctv, col = collist)
hm.celltype <- Heatmap(matrix = hm2,
                       show_column_dend = FALSE,
                       show_row_dend = FALSE,
                       name = hm.name,
                       rect_gp = gpar(col = "black", lwd = 1))
  
# heatmap -- composite
fname <- paste0("hm-markers-composite_sample-celltype_marker-", 
                celltype.variable,"_", proj.str,".jpg")
jpeg(fname, width = 8, height = 7, units = "in", res = 400)
print(hm.sample.celltype + hm.celltype); dev.off()

