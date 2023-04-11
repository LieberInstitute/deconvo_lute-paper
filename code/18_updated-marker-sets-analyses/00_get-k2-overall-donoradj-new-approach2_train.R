
libv <- c("here", "dplyr", "lute", "biomaRt", "SummarizedExperiment", 
          "SingleCellExperiment", "scuttle")
sapply(libv, library, character.only = TRUE)

# main run parameters
markers.per.group.discover <- 500
markers.per.group.final <- 40
celltype.variable.original <- "cellType_broad_hc"
celltype.variable.new <- "k2"
group.id.variable <- "Sample"
assay.name.markers <- "logcounts"
typemarker.algorithm <- "meanratios"
min.group.overlap.rate <- 0.8
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
# get all marker candidates
markers.by.group <- markers_by_group(sce, 
                                     group.variable = group.id.variable, 
                                     celltype.variable = celltype.variable.new, 
                                     assay.name = assay.name.markers, 
                                     markers.per.type = markers.per.group.discover, 
                                     typemarker.algorithm = typemarker.algorithm,
                                     return.type = "list",
                                     verbose = TRUE)
# get marker overlap info
markers.filtered <- filter_group_markers(markers.by.group,
                                         minimum.group.overlap.rate = 
                                           min.group.overlap.rate)
# save overlap markers info
markers.filtered.path <- file.path("deconvo_method-paper/outputs/", 
                                   markers.filtered.name)
save(markers.filtered, file = markers.filtered.path)

# filter marker se data

# get slide adj counts

# get log counts

# get mean ratios results

# get signature matrix






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
utypev <- unique(sce[[celltype.variable.new]])
sce.ds <- do.call(cbind, lapply(utypev, function(ti){
  message("downsampling for type ", ti, "...")
  # filter sce
  filter.sce <- sce.filter[[celltype.variable.new]]==ti
  scef <- sce.filter[,filter.sce]
  # get filtered data
  batchv <- scef[[batch.variable.name]]
  mexpr <- assays(scef)[[assay.name.adj]]
  # downsample
  mexpr.ds <- downsampleBatches(mexpr, batch = batchv)
  assays(scef)[[assay.name.adj]] <- mexpr.ds 
  scef
}))
message("Getting top markers on downsampled log slide-adjusted counts..")
sce.ds <- logNormCounts(sce.ds, assay.type = assay.name.adj)
save(sce.ds, file = 
       here("deconvo_method-paper/outputs/sce-all_slide-adj-ds_train.rda"))

#--------------------------------
# get markers from slide-adj data
#--------------------------------
# get markers
markers.results.adj.1k <- meanratiosParam(sce.ds, 
                                          assay.name.adj, 
                                          celltype.variable.new, 
                                          markers.per.type = 500) %>%
  typemarkers()
markers.results.adj.40 <- meanratiosParam(sce.ds, 
                                          assay.name.adj, 
                                          celltype.variable.new, 
                                          markers.per.type = 20) %>%
  typemarkers()
# make new se
filter.sce <- rownames(sce.ds) %in% markers.results.adj$allmarkers
sce.ds.filter <- sce.ds[filter.sce,]
se.filter <- SummarizedExperiment(assays = assays(sce.ds.filter)[["counts"]],
                                  colData = colData(sce.ds.filter),
                                  rowData = rowData(sce.ds.filter))
metadata(se.filter) <- list(top.500.markers.results = markers.results.adj.1k,
                            top.40.markers.results = markers.results.adj.40)
# save
se.filter.all.adj.name <- "se-adj_marker-filter-all_train.rda"
save(se.filter, file = here(save.path, se.filter.all.adj.name))

#-----------------------------------------------------
# get overlaps with k2 concordant, overlapping markers
#-----------------------------------------------------



