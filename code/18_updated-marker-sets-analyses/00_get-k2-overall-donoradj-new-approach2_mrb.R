# k2 markers in mrb dataset, including with donor adj counts/logcounts

libv <- c("here", "dplyr", "lute", "biomaRt", "SummarizedExperiment", 
          "SingleCellExperiment", "scuttle", "sva")
sapply(libv, library, character.only = TRUE)

# main run parameters
markers.per.group.discover <- 500
markers.per.group.final <- 40
celltype.variable.original <- "cellType"
celltype.variable.new <- "k2"
group.id.variable <- "donor"
assay.name.markers <- "counts"
typemarker.algorithm <- "meanratios"
min.group.overlap.rate <- 0.8
markers.bygroup.name <- "group-markers-k2_mrb.rda"
markers.filtered.name <- "overlap-markers-filtered-k2_mrb.rda"

# load sce
sce.mrb.name <- "sce-mrb_dlpfc.rda"
sce.mrb.path <- here("deconvo_method-paper", "outputs", "09_manuscript", sce.mrb.name)
sce <- get(load(sce.mrb.path))

# get celltype variable
celltype.variable.vector <- ifelse(grepl("Excit|Inhib", sce[["cellType"]]),"neuron",
                                   ifelse(sce[["cellType"]] %in% 
                                            c("Oligo", "OPC", "Micro", "Astro"), 
                                          "glial", "other"))
sce[[celltype.variable.new]] <- celltype.variable.vector
# remove "other"
filter.sce <- !sce[[celltype.variable.new]] == "other"
sce <- sce[,filter.sce]

sce <- scuttle::logNormCounts(sce)

# get all marker candidates
markers.by.group <- markers_by_group(sce, 
                                     group.variable = group.id.variable,
                                     celltype.variable = celltype.variable.new, 
                                     assay.name = "logcounts", 
                                     markers.per.type = markers.per.group.discover, 
                                     typemarker.algorithm = typemarker.algorithm,
                                     return.type = "list", verbose = TRUE)
# save
markers.bygroup.path <- file.path("deconvo_method-paper/outputs/", 
                                  markers.bygroup.name)
save(markers.by.group, file = markers.bygroup.path)

# marker overlaps/concordance filter
markers.filtered <- filter_group_markers(markers.by.group,
                                         minimum.group.overlap.rate = 0.5)
# save overlap markers info
markers.filtered.path <- file.path("deconvo_method-paper/outputs/", 
                                   markers.filtered.name)
save(markers.filtered, file = markers.filtered.path)

#---------------------
# get slide-adj counts
#---------------------
dim(sce) # [1] 36601 77604
scale.thresh <- -0.9 # scale filter
assay.name.adj <- "counts_adjusted"
celltype.variable.new <- "k2"
batch.variable.name <- "donor"

# get downsampled cell type amounts
cd <- colData(sce)
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
dim(sce.filter) # [1] 33538 11165

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
# Found 7428 genes with uniform expression within a single batch (all zeros); 
# these will not be adjusted for batch

# save sce adj
save(sce.filter, file = 
       here("deconvo_method-paper/outputs/sce-all_slide-adj_mrb.rda"))

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
message("Getting top markers on downsampled log slide-adjusted counts..")
sce.ds <- scuttle::logNormCounts(sce.ds, assay.type = assay.name.adj)
save(sce.ds, file = 
       here("deconvo_method-paper/outputs/sce-all_slide-adj-ds_mrb.rda"))

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
se.filter.all.adj.name <- "se-adj_marker-filter-all_mrb.rda"
save(se.filter, file = 
       here("deconvo_method-paper", "outputs", se.filter.all.adj.name))

#-----------------------------------------------------
# get overlaps with k2 concordant, overlapping markers
#-----------------------------------------------------
markers.concordant <- get(load(here("deconvo_method-paper", 
                                    "outputs", "overlap-markers-filtered-k2_train.rda")))
length(intersect(rownames(se.filter), markers.concordant[[1]]))
# [1] 197

