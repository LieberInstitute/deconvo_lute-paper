# make formatted heatmap of new k2 markers.
# use the concordant, overlapping marker subset from training data.
# 

library(ComplexHeatmap)

outputs.path <- here("deconvo_method-paper", "outputs", "18_updated-marker-sets-analyses")

se.concordant.path <- here(outputs.path, 
                           "se-adj_overlap-concordant_marker-filter-all_train.rda")
se.train.concordant <- get(load(se.concordant.path))

#-------------
# prep se data
#-------------
# filter se
se <- se.train.concordant
matrix.counts <- assays(se)[[assay.name]]
matrix.counts[matrix.counts < 0] <- 0
filter.zeros <- colSums(matrix.counts) == 0
dim(matrix.counts)
matrix.counts <- matrix.counts[,!filter.zeros]
dim(matrix.counts)
filter.se <- colnames(se) %in% colnames(matrix.counts)
se.filter <- se[,filter.se]
se.filter <- scuttle::logNormCounts(se.filter, assay.type = "counts_adj")
# convert to sce
sce <- SingleCellExperiment(se.filter)
colData(sce) <- colData(se.filter)
assays(sce) <- list(counts = as.matrix(assays(se.filter)[["counts_adj"]]),
                    logcounts = as.matrix(assays(se.filter)[["logcounts"]]))

#-------------------
# get top 40 markers
#-------------------
top.markers <- meanratiosParam(sce, assay.name = "counts",
                               celltype.variable = "k2", markers.per.type = 20,
                               return.info = T) %>% typemarkers()

filter.sce <- rownames(sce) %in% top.markers
sce.filter <- sce[filter.sce,]
dim(sce.filter)

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
sce.hm <- sce.filter
# summarize logcounts
celltype.vector <- sce.hm[[celltype.variable]]
group.vector <- sce.hm[[group.variable]]
new.group.label.vector <- paste0(group.vector, ";", celltype.vector)
sce.means <- aggregateAcrossCells(x = sce.hm, ids = new.group.label.vector,
                                  use.assay.type = "logcounts")
# scale expression
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
sce.means2 <- aggregateAcrossCells(x = sce.hm, ids = celltype.vector,
                                   use.assay.type = "logcounts")
me2 <- assays(sce.means2)[["logcounts"]]
cnv <- colnames(me2)
marker.vector <- unlist(apply(me2, 1, function(ri){
  cnv[which(ri==max(ri))]}))
marker.vector <- as.character(marker.vector)
leftanno <- rowAnnotation(celltype = marker.vector,
                          col = list(celltype = c("neuron" = "blue",
                                                  "glial" = "yellow")),
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
fname <- paste0("hm-markers_overlapping-concordant_sample-celltype_marker-", 
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
fname <- paste0("hm-markers-composite_concordant-overlapping_sample-celltype_marker-", 
                celltype.variable,"_", proj.str,".jpg")
jpeg(fname, width = 8, height = 7, units = "in", res = 400)
print(hm.sample.celltype + hm.celltype); dev.off()

#-------------------------
# heatmap with old markers
#-------------------------
