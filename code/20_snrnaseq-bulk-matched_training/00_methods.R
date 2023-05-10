# dependencies
libv <- c("here", "lute", "dplyr", "ggplot2", "gridExtra", "SingleCellExperiment", 
          "SummarizedExperiment", "scran", "DeconvoBuddies", "UpSetR", "DelayedArray")
sapply(libv, library, character.only = TRUE)

# save path
save.path <- here("deconvo_method-paper", "outputs", "20_snrnaseq-bulk-matched_training")

#--------------------------
# get markers within donors
#--------------------------
# we selected the top markers within each donor
num.markers <- 500

# load sn expression
sce.path <- file.path("DLPFC_snRNAseq/processed-data/sce","sce_DLPFC.Rdata")
sce <- get(load(sce.path))
# filter training samples

# for each donor, get the top markers
cell.type.variable <- "cellType_broad_hc"
sample.id.vector <- unique(sce$Sample)
k2.types.labels <- list("neuron" = c("Excit", "Inhib"), 
                        "glial" = c("OPC", "Oligo", "Astro", "Micro"))
k3.types.labels <- list("neuron" = c("Excit", "Inhib"), "oligo" = c("Oligo"),
                        "glial_non_oligo" = c("Astro", "Micro"))
k4.types.labels <- list("neuron" = c("Excit", "Inhib"), "oligo" = c("Oligo"),
                        "astro" = c("Astro"), "micro" = c("Micro"))
sce.filter <- sce[,sce[[cell.type.variable]] %in% unlist(k2.types.labels)]
# apply labels


list.markers.sn <- lapply(sample.id.vector, function(sample.id){
  # k2 
  sce.filter <- sce[[cell.type.variable]] %in% unlist(k2.types.labels)
  sce.filter.k2 <- sce[, ]
  markers.k2 <- lute(
    sce = sce.filter.k2,
    markers.per.type = num.markers, 
       celltype.variable = "k2", 
       deconvolution.algorithm = NULL)
  # k3
  # k4
})

# select markers within every matched sample

#-----------------
# matched datasets
#-----------------
#sce.prepared.path <- file.path("deconvo_method-paper", "outputs", "09_manuscript",
#                               "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda")
#lsce <- get(load(sce.prepared.path))

# get sce for varying k values
#sce.k2 <- lsce$k2
#sce.k3 <- lsce$k3
#sce.k4 <- lsce$k4

# bulk rnaseq
# paths
rse.bulk.filename <- "rse_gene.Rdata"
rse.bulk.filepath <- here("Human_DLPFC_Deconvolution", "processed-data", 
                          "01_SPEAQeasy", "round2_v40_2022-07-06", "rse", 
                          rse.bulk.filename)
rse <- get(load(rse.bulk.filepath))
rse$Sample <- paste0(rse$BrNum, "_", tolower(rse$location))

# assemble matched experiments
sample.id.vector.sn <- unique(sce.k2$Sample)
sample.id.vector.bulk <- unique(rse$Sample)
decon.k2 <- lapply(sample.id.vector.sn, function(sample.id){
  sce.sample <- sce[,sce$Sample==sample.id]
  rse.sample <- rse[,rse$Sample == sample.id]
  if(ncol(rse.sample)>0){
    # get markers 
    list(sce = sce.sample, y.series = rse.sample)
  }
})
decon.k3 <- lapply()
decon.k4 <- lapply()

# get deconvolution results for matched samples
# without sfactor transformation
s <- c(1, 1)

# get sce for varying k values
#sce.k2 <- lsce$k2
#sce.k3 <- lsce$k3
#sce.k4 <- lsce$k4
# with some sfactor transformation
s <- c(3, 10)
