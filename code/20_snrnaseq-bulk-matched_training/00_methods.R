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
cell.type.variable <- "cellType_broad_hc"

# load sn expression
sce.path <- file.path("DLPFC_snRNAseq/processed-data/sce","sce_DLPFC.Rdata")
sce <- get(load(sce.path))
# filter training samples
donor.variable <- "BrNum"
donor.id.train <- c("Br2720", "Br6471", "Br8492", "Br2743", "Br3942", "Br6423", "Br8325")
filter.train <- sce[[donor.variable]] %in% donor.id.train
sce <- sce[,filter.train]

# for each donor, get the top markers
sample.id.vector <- unique(sce$Sample)
k2.types.labels <- list("neuron" = c("Excit", "Inhib"), 
                        "glial" = c("OPC", "Oligo", "Astro", "Micro"))
k3.types.labels <- list("neuron" = c("Excit", "Inhib"), "oligo" = c("Oligo"),
                        "glial_non_oligo" = c("Astro", "Micro"))
k4.types.labels <- list("neuron" = c("Excit", "Inhib"), "oligo" = c("Oligo"),
                        "astro" = c("Astro"), "micro" = c("Micro"))
sce <- sce[,sce[[cell.type.variable]] %in% unlist(k2.types.labels)]
# apply labels
sce[["k2"]] <- ifelse(sce[[cell.type.variable]] %in% k2.types.labels[["neuron"]], "neuron", "glial")
sce[["k3"]] <- ifelse(sce[[cell.type.variable]] %in% k3.types.labels[["neuron"]], "neuron", 
                      ifelse(sce[[cell.type.variable]] %in% k3.types.labels[["oligo"]], 
                             "oligo", "glial_non_oligo"))
sce[["k4"]] <- ifelse(sce[[cell.type.variable]] %in% k4.types.labels[["neuron"]], "neuron", 
                      ifelse(sce[[cell.type.variable]] %in% k4.types.labels[["oligo"]], "oligo", 
                             ifelse(sce[[cell.type.variable]] %in% k4.types.labels[["astro"]], 
                                    "astro", "micro")))
# save
sce.prep.path <- './deconvo_method-paper/outputs/sce_prep_05-10-23_train.rda'
save(sce, file = sce.prep.path)

num.markers <- 100
list.markers.sn <- lapply(sample.id.vector, function(sample.id){
  sce.filter <- sce[,sce$Sample==sample.id]
  # k2
  markers.k2 <- lute(
    sce = sce.filter,
    markers.per.type = num.markers, 
       celltype.variable = "k2", 
       deconvolution.algorithm = NULL)
  # k3
  markers.k3 <- lute(
    sce = sce.filter,
    markers.per.type = num.markers, 
    celltype.variable = "k3", 
    deconvolution.algorithm = NULL)
  # k4
  markers.k4 <- lute(
    sce = sce.filter,
    markers.per.type = num.markers, 
    celltype.variable = "k4", 
    deconvolution.algorithm = NULL)
  return(list(k2 = markers.k2, k3 = markers.k3, k4 = markers.k4))
})
names(list.markers.sn) <- sample.id.vector
list.markers.path <- "./deconvo_method-paper/outputs/list-markers-by-donor-k2-3-4_train.rda"
save(list.markers.sn, file = list.markers.path)

# subset sce on only markers
list.markers.path <- "./deconvo_method-paper/outputs/list-markers-by-donor-k2-3-4_train.rda"
list.markers.sn <- get(load(list.markers.path))
all.markers.vector <- lapply(list.markers.sn, function(list.iter){
  c(list.iter[[1]][[1]], list.iter[[2]][[1]], list.iter[[3]][[1]])
}) %>% unlist() %>% unique()
sce <- sce[rownames(sce) %in% all.markers.vector,]

#-----------------------------------
# do deconvolution with matched data
#-----------------------------------
# bulk rnaseq
# paths
rse.path <- "Human_DLPFC_Deconvolution/processed-data/rse/rse_gene.Rdata"
rse <- get(load(rse.path))
rse$Sample <- paste0(rse$BrNum, "_", tolower(rse$location))
rse <- rse[!duplicated(rowData(rse)$Symbol),]
rownames(rse) <- rowData(rse)$Symbol

# get experiment results
decon.all <- lapply(sample.id.vector.sn, function(sample.id){
  sce.sample <- sce[,sce$Sample == sample.id]
  rse.sample <- rse[,rse$Sample == sample.id]
  if(ncol(rse.sample) > 0){
    for(markers in list.markers.sn[[sample.id]]){
      sce.markers <- sce.sample[rownames(sce.sample) %in% markers[[1]],]
      rse.markers <- rse.sample[rownames(rse.sample) %in% markers[[1]],]
      # without s
      decon.no.scale <- lute()
      # with s
      decon.with.scale <- lute()
    }
  }
})


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
