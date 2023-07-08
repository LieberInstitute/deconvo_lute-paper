
libv <- c("here", "dplyr", "lute", "biomaRt", "SummarizedExperiment", 
          "SingleCellExperiment", "scuttle")
sapply(libv, library, character.only = TRUE)

# main run parameters
markers.per.group.discover <- 500
markers.per.group.final <- 40
celltype.variable.original <- "cellType_broad_hc"
celltype.variable.new <- "k2"
group.id.variable <- "Sample"
assay.name.markers <- "counts"
typemarker.algorithm <- "meanratios"
min.group.overlap.rate <- 0.5
markers.bygroup.name <- "group-markers-k2_train.rda"
markers.filtered.name <- "overlap-markers-filtered-k2_train.rda"

# get sce data
sce.name <- "sce_DLPFC.Rdata"
sce.path <- file.path("DLPFC_snRNAseq/processed-data/sce", sce.name)
sce.prepared.path <- file.path("deconvo_method-paper", "outputs", 
                               "15_sample-wise-signature-matrix-simulations",
                               "sce-prepared_dlpfc-ro1-train.rda")
sce <- get(load(sce.path))
# subset training donors
donor.variable <- "BrNum"
donor.id.train <- c("Br2720", "Br6471", "Br8492", "Br2743", "Br3942", "Br6423", "Br8325")
filter.train <- sce[[donor.variable]] %in% donor.id.train
sce <- sce[,filter.train]
# get k2 labels and filter
cell.type.vector <- sce[[celltype.variable.original]]
k2.cell.type.vector <- ifelse(cell.type.vector %in% c("Excit", "Inhib"), "neuron",
                              ifelse(cell.type.vector %in% c("Oligo", "OPC", "Astro", "Micro"), 
                                     "glial", "other"))
sce[[celltype.variable.new]] <- k2.cell.type.vector
filter.other <- !sce[[celltype.variable.new]]=="other"
sce <- sce[,filter.other]
sce <- scuttle::logNormCounts(sce)

#scef <- sce[,sce[["Sample"]]=="Br2720_post"]
#markers <- meanratiosParam(scef, "counts", "k2", 500, TRUE) %>% typemarkers()

# get all marker candidates
# omit 2720
sce <- sce[,!sce[[group.id.variable]] == "Br2720_post"]
markers.by.group <- markers_by_group(sce, 
                                     group.variable = group.id.variable, 
                                     celltype.variable = celltype.variable.new, 
                                     assay.name = "logcounts", 
                                     markers.per.type = markers.per.group.discover, 
                                     typemarker.algorithm = typemarker.algorithm,
                                     return.type = "list",
                                     verbose = TRUE)

# save
markers.bygroup.path <- file.path("deconvo_method-paper/outputs/", 
                                   markers.bygroup.name)
save(markers.by.group, file = markers.bygroup.path)

#-------------------------------
# visualize concordance overlaps
#-------------------------------

# upset plot of markers by participant
library(UpSetR)
unique.types <- c("glial", "neuron")
unique.donors <- names(markers.by.group)
list.concord <- lapply(unique.types, function(type){
  lapply(unique.donors, function(donor){
    table <- markers.by.group[[donor]]
    tf <- table[table$cellType.target==type,]
    tf$gene
  })
}) %>% unlist(recursive = F)
names(list.concord) <- paste0(unique.donors, ";", 
                              rep(unique.types, each = length(unique.donors)))
upset(fromList(list.concord), nsets = 20)

# barplots
# get plot data
glial.all <- list.concord[grepl("glial", names(list.concord))] %>% unlist()
neuron.all <- list.concord[grepl("neuron", names(list.concord))] %>% unlist()
nonconcord.all <- intersect(glial.all, neuron.all)
length(nonconcord.all) # 135
concord.glial <- glial.all[!glial.all %in% nonconcord.all]
length(concord.glial) # 136
concord.neuron <- neuron.all[!neuron.all %in% nonconcord.all]
length(concord.neuron) # 5297
# percent concordant
length(concord.neuron)/length(neuron.all) # 0.9345448
length(concord.glial)/length(glial.all) # 0.4788732
# mean overlapping
glial.overlap <- list.concord[grepl("glial", names(list.concord))] %>% 
  unlist() %>% table() %>% as.data.frame()
glial.overlap <- glial.overlap[,2]
glial.mean.overlap <- mean(glial.overlap) %>% round(digits = 2) # 1.083969
glial.sd.overlap <- sd(glial.overlap) %>% round(digits = 2)
glial.median.overlap <- median(glial.overlap) %>% round(digits = 2)
neuron.overlap <- list.concord[grepl("neuron", names(list.concord))] %>% 
  unlist() %>% table() %>% as.data.frame()
neuron.overlap <- neuron.overlap[,2]
neuron.mean.overlap <- mean(neuron.overlap) %>% round(digits = 2)
neuron.sd.overlap <- sd(neuron.overlap) %>% round(digits = 2)
neuron.median.overlap <- median(neuron.overlap) %>% round(digits = 2)
data.frame(type = c("neuron", "glial"),
           mean = c(neuron.mean.overlap, glial.mean.overlap),
           sd = c(neuron.sd.overlap, glial.sd.overlap),
           median = c(neuron.median.overlap, glial.median.overlap))
# make ggplot
plot.data <- data.frame(marker.type = c("glial", "neuron", "both"),
                        count = c(length(concord.glial), 
                                  length(concord.neuron), 
                                  length(nonconcord.all)))
plot.data[,2] <- as.numeric(plot.data[,2])
barplot.markers <- ggplot(plot.data, aes(x = marker.type, y = count)) + 
  geom_bar(stat = "identity") + geom_label(aes(y=count,label=count)) + 
  scale_y_log10() + ylab("Marker count (log10 scaled)") + theme_bw()
# save plot
jpeg("barplot_marker-counts-by-type_train.jpg", width = 3, height = 3, units = "in", res = 400)
barplot.markers; dev.off()


# get marker overlap info
# get marker table list
index.vector <- seq(length(markers.by.group))
list.marker.tables <- lapply(index.vector, function(index){
  marker.table <- markers.by.group[[index]]$result.info
  marker.table$group.id <- names(markers.by.group)[index]
  return(marker.table)
})
markers.filtered <- filter_group_markers(list.marker.tables,
                                         minimum.group.overlap.rate = 
                                           min.group.overlap.rate)
# Found 8085 total markers.
# Getting marker overlap info...
# Found 4400 concordant markers by type.
# Found 99 concordant markers with overlap rate at least 0.5

# save overlap markers info
markers.filtered.path <- file.path("deconvo_method-paper/outputs/",markers.filtered.name)
save(markers.filtered, file = markers.filtered.path)

#---------------------
# get slide-adj counts
#---------------------
libv <- c("here", "dplyr", "lute", "biomaRt", "SummarizedExperiment", 
          "SingleCellExperiment", "scuttle", "sva")
sapply(libv, library, character.only = TRUE)

# get sce data
sce.name <- "sce_DLPFC.Rdata"
sce.path <- file.path("DLPFC_snRNAseq/processed-data/sce", sce.name)
sce <- get(load(sce.path))

filter.sce <- rownames(sce) %in% markers.filtered$marker.list.final

dim(sce) # [1] 36601 77604

scale.thresh <- -0.9 # scale filter
assay.name.adj <- "counts_adjusted"
celltype.variable.original <- "cellType_broad_hc"
celltype.variable.new <- "k2"
batch.variable.name <- "Sample"

# subset training donors
donor.variable <- "BrNum"
donor.id.train <- c("Br2720", "Br6471", "Br8492", "Br2743", "Br3942", "Br6423", "Br8325")
filter.train <- sce[[donor.variable]] %in% donor.id.train
sce <- sce[,filter.train]
# get k2 labels and filter
cell.type.vector <- sce[[celltype.variable.original]]
k2.cell.type.vector <- ifelse(cell.type.vector %in% c("Excit", "Inhib"), "neuron",
                              ifelse(cell.type.vector %in% c("Oligo", "OPC", "Astro", "Micro"), 
                                     "glial", "other"))
sce[[celltype.variable.new]] <- k2.cell.type.vector
filter.other <- !sce[[celltype.variable.new]]=="other"
sce <- sce[,filter.other]
sce <- scuttle::logNormCounts(sce)
cd <- colData(sce)

# get downsampled cell type amounts
dfm <- cd[,celltype.variable.new] %>% table() %>% as.data.frame()
type.filt <- as.character(dfm[dfm[,2]==min(dfm[,2]),1])
cdf <- cd[cd[,celltype.variable.new] == type.filt,]
dft <- table(cdf[,batch.variable.name], cdf[,celltype.variable.new]) %>% as.data.frame()
scalev <- scale(dft$Freq)[,1]
which.filt <- which(scalev <= scale.thresh)
sample.filt <- as.character(dft[which.filt,1])
filter.cd <- cd[,batch.variable.name] %in% sample.filt
cells.remove <- cd[filter.cd,] %>% rownames()
sce.filter <- sce[,!colnames(sce) %in% cells.remove]
md.rm <- list(type.filt = type.filt,
              scale.thresh = scale.thresh,
              num.cells = length(cells.remove),
              batch.id = sample.filt)
message("removed ", length(cells.remove),
        " cells from ", length(sample.filt),
        " batches for cell type: '", type.filt, "'.")
# removed 4997 cells from 2 batches for cell type: 'glial'.
dim(sce.filter)

message("doing combat adj...")
mexpr <- assays(sce.filter)[["counts"]]
cnv <- colnames(sce.filter)
pheno <- data.frame(donor = sce.filter[[batch.variable.name]], 
                    type = sce.filter[[celltype.variable.new]])
mod <- model.matrix(~type, data = pheno)
mi.adj <- ComBat(dat = mexpr, batch = pheno$donor, mod = mod)
message("converting negative values...")
mi.adj[mi.adj < 0] <- 0 # convert negative values
assays(sce.filter)[[assay.name.adj]] <- mi.adj
# note: 10027 / 36601 genes have uniform variances, so aren't adjusted
# save sce adj
save(sce.filter, file = 
       here("deconvo_method-paper/outputs/sce-all_slide-adj_train.rda"))

message("performing downsampling...")
# downsample -- by donor, within celltypes
utypev <- unique(sce.filter[[celltype.variable.new]])
sce.ds <- do.call(cbind, lapply(utypev, function(ti){
  message("downsampling for type ", ti, "...")
  # filter sce
  filter.sce <- sce.filter[[celltype.variable.new]]==ti
  scef <- sce.filter[,filter.sce]
  # get filtered data
  batchv <- scef[[batch.variable.name]]
  mexpr <- assays(scef)[[assay.name.adj]]
  # downsample
  mexpr.ds <- scuttle::downsampleBatches(mexpr, batch = batchv)
  assays(scef)[[assay.name.adj]] <- mexpr.ds 
  scef
}))
sce.ds <- scuttle::logNormCounts(sce.ds, assay.type = assay.name.adj)
#save(sce.ds, file = 
#       here("deconvo_method-paper/outputs/sce-all_slide-adj-ds_train.rda"))

# filter overlapping markers
filter.sce.ds <- rownames(sce.ds) %in% markers.filter$marker.list.final
sce.ds.filtered <- sce.ds[filter.sce.ds,]
# make new se
matrix.counts.adj <- assays(sce.ds.filtered)[["counts_adjusted"]] %>% as.matrix()
assays.list <- list(counts_adj = matrix.counts.adj)
se.filter <- SummarizedExperiment(assays = assays.list,
                                  colData = colData(sce.ds.filtered),
                                  rowData = rowData(sce.ds.filtered))
metadata(se.filter) <- list(markers.by.group = markers.by.group,
                            markers.overlap.info = markers.filter)
# save new se
se.filter.all.adj.name <- "se-adj_overlap-concordant_marker-filter-all_train.rda"
save(se.filter, file = here("deconvo_method-paper", "outputs", se.filter.all.adj.name))

#--------------------------------
# get markers from slide-adj data
#--------------------------------
# get markers
markers.results.adj.1k <- meanratiosParam(sce.ds, 
                                          "logcounts", 
                                          celltype.variable.new, 
                                          markers.per.type = 500,
                                          return.info = TRUE) %>%
  typemarkers()
markers.results.adj.40 <- meanratiosParam(sce.ds, 
                                          "logcounts", 
                                          celltype.variable.new, 
                                          markers.per.type = 20,
                                          return.info = TRUE) %>%
  typemarkers()
# make new se
filter.sce <- rownames(sce.ds) %in% markers.results.adj.1k$markers
sce.ds.filter <- sce.ds[filter.sce,]
matrix.counts <- assays(sce.ds.filter)[["counts"]] %>% as.matrix()
matrix.logcounts <- assays(sce.ds.filter)[["logcounts"]] %>% as.matrix()
matrix.counts.adj <- assays(sce.ds.filter)[["counts_adjusted"]] %>% as.matrix()
#assays.list <- list(counts = matrix.counts,
#                    logcounts_adj = matrix.logcounts,
#                    counts_adj = matrix.counts.adj)
assays.list <- list(counts_adj = matrix.counts.adj)
se.filter <- SummarizedExperiment(assays = assays.list,
                                  colData = colData(sce.ds.filter),
                                  rowData = rowData(sce.ds.filter))
metadata(se.filter) <- list(top.500.markers.results = markers.results.adj.1k,
                            top.40.markers.results = markers.results.adj.40)
# save
se.filter.all.adj.name <- "se-adj_marker-filter-all_train.rda"
save(se.filter, file = here("deconvo_method-paper", "outputs", se.filter.all.adj.name))

#-----------------------------------------------------
# get overlaps with k2 concordant, overlapping markers
#-----------------------------------------------------
markers.concordant <- get(load(here("deconvo_method-paper", 
                                    "outputs", "overlap-markers-filtered-k2_train.rda")))
length(intersect(rownames(se.filter), markers.concordant[[1]]))
# [1] 197

