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
# load list markers
list.markers.sn <- get(load("./deconvo_method-paper/outputs/list-markers-by-donor-k2-3-4_train.rda"))
# load sce
sce <- get(load('./deconvo_method-paper/outputs/sce_prep_05-10-23_train.rda'))
# bulk rnaseq
# paths
rse.path <- "Human_DLPFC_Deconvolution/processed-data/rse/rse_gene.Rdata"
rse <- get(load(rse.path))
# rse$Sample <- paste0(rse$BrNum, "_", ifelse(rse))
rse <- rse[!duplicated(rowData(rse)$Symbol),]
rownames(rse) <- rowData(rse)$Symbol

# get experiment results
s.vector <- c("glial" = 3, "glial_non_oligo" = 3, 
              "neuron" = 10, "astro" = 4, "oligo" = 3, "micro" = 3)
sample.id.vector.sn <- unique(sce$Sample)
result.list <- lapply(sample.id.vector.sn, function(sample.id){
  sce.sample <- sce[,sce$Sample == sample.id]
  rse.sample <- rse[,rse$Sample == sample.id]
  if(ncol(rse.sample) > 0){
    sample.markers <- list.markers.sn[[sample.id]]
    marker.index.vector <- seq(length(sample.markers))
    results.bymarker.list <- lapply(marker.index.vector, function(marker.index){
      # filter and match marker order
      marker.type <- names(sample.markers)[marker.index]
      markers <- sample.markers[[marker.index]]
      markers.vector <- markers$typemarker.results
      sce.markers <- sce.sample[rownames(sce.sample) %in% markers.vector,]
      rse.markers <- rse.sample[rownames(rse.sample) %in% markers.vector,]
      shared.markers <- intersect(rownames(sce.markers), rownames(rse.markers))
      sce.markers <- sce.sample[rownames(sce.sample) %in% shared.markers,]
      rse.markers <- rse.sample[rownames(rse.sample) %in% shared.markers,]
      rse.markers <- rse.markers[order(match(rownames(rse.markers), shared.markers)),]
      identical(rownames(rse.markers), rownames(sce.markers))
      
      # make y expression
      y.expression <- assays(rse.markers)[[1]]
      colnames(y.expression) <- paste0(rse.markers$Sample, "_", 
                                       rse.markers$library_type, "_", 
                                       rse.markers$library_prep)
      
      # without s
      decon.no.scale <- lute(sce = sce.markers, y = y.expression,
                             celltype.variable = marker.type,
                             typemarker.algorithm = NULL)$deconvolution.results %>%
        as.data.frame()
      # with s
      decon.with.scale <- lute(sce = sce.markers, y = y.expression,
                               celltype.variable = marker.type, s = s.vector,
                               typemarker.algorithm = NULL)$deconvolution.results %>%
        as.data.frame()
      decon.no.scale$scale <- FALSE
      decon.with.scale$scale <- TRUE
      df.result <- rbind(decon.no.scale, decon.with.scale) %>% as.data.frame()
      df.result$sample_label <- rownames(df.result)
      df.result$marker.type <- marker.type
      return(df.result)
    })
    # results.bymarker.table <- do.call(rbind, results.bymarker.list) %>% as.data.frame()
    return(results.bymarker.list)
  }
})
result.table.k2 <- do.call(rbind, 
                           lapply(result.list, 
                                  function(list.iter){list.iter[[1]]})) %>% as.data.frame()
result.table.k3 <- do.call(rbind, 
                           lapply(result.list, 
                                  function(list.iter){list.iter[[2]]})) %>% as.data.frame()
result.table.k4 <- do.call(rbind, 
                           lapply(result.list, 
                                  function(list.iter){list.iter[[3]]})) %>% as.data.frame()

deconvo.results.list.k234 <- list(k2 = result.table.k2,
                                  k3 = result.table.k3,
                                  k4 = result.table.k4)

result.table.name <- "deconvo-results-table_bulk-matched-sn-bydonor_train.rda"
result.table.path <- paste0("./deconvo_method-paper/outputs/", result.table.name)
save(deconvo.results.list.k234, file = result.table.path)

#---------------
# map image data
#---------------
# load data
# deconvo results table
result.table.list <- get(load("./deconvo_method-paper/outputs/20_snrnaseq-bulk-matched_training/deconvo-results-table_bulk-matched-sn-bydonor_train.rda"))
# image reference
halo.path <- "Human_DLPFC_Deconvolution/processed-data/03_HALO/halo_all.Rdata"
halo.all <- get(load(halo.path))

# cell type proportions by sample
halo.sample.prop <- do.call(rbind, sapply(unique(halo.all$Sample), function(sample.id){
  table(halo.all[halo.all$Sample==sample.id,]$cell_type) %>% prop.table()
})) %>% as.data.frame()
halo.sample.prop$sample.id <- unique(halo.all$Sample)

#-----------------
# merge k2 results
#-----------------
result.filter <- result.table.list["k2"][[1]]
result.filter$halo.neuron <- result.filter$halo.glial <- 0
for(ii in seq(nrow(result.filter))){
  sample.id <- result.filter$sample_label[ii]
  sample.id.format <- paste0(unlist(strsplit(sample.id, "_"))[1:2], collapse = "_")
  which.halo.id <- halo.sample.prop$sample.id==sample.id.format
  halo.filter <- halo.sample.prop[which.halo.id,]
  result.filter$halo.neuron[ii] <- halo.filter$Excit+halo.filter$Inhib
  result.filter$halo.glial[ii] <- halo.filter$Oligo+halo.filter$Astro+halo.filter$Micro
}
# get errors
result.filter$error.neuron <- result.filter$halo.neuron-result.filter$neuron
result.filter$abs.error.neuron <- abs(result.filter$error.neuron)
result.filter$error.glial <- result.filter$halo.glial-result.filter$glial
result.filter$abs.error.glial <- abs(result.filter$glial)

#-----------------
# merge k3 results
#-----------------
result.filter <- result.table.list["k3"][[1]]
result.filter$halo.neuron <- result.filter$halo.oligo <- result.filter$halo.non.oligo.glial <- 0
for(ii in seq(nrow(result.filter))){
  sample.id <- result.filter$sample_label[ii]
  sample.id.format <- paste0(unlist(strsplit(sample.id, "_"))[1:2], collapse = "_")
  which.halo.id <- halo.sample.prop$sample.id==sample.id.format
  halo.filter <- halo.sample.prop[which.halo.id,]
  result.filter$halo.neuron[ii] <- halo.filter$Excit+halo.filter$Inhib
  result.filter$halo.oligo[ii] <- halo.filter$Oligo
  result.filter$halo.non.oligo.glial[ii] <- halo.filter$Astro+halo.filter$Micro
}
# get errors
result.filter$error.neuron <- result.filter$halo.neuron-result.filter$neuron
result.filter$abs.error.neuron <- abs(result.filter$error.neuron)
result.filter$error.oligo <- result.filter$halo.oligo-result.filter$oligo
result.filter$abs.error.oligo <- abs(result.filter$oligo)
result.filter$error.non.oligo.glial <- result.filter$halo.non.oligo.glial-result.filter$glial_non_oligo
result.filter$abs.error.non.oligo.glial <- abs(result.filter$error.non.oligo.glial)

#-----------------
# merge k4 results
#-----------------
result.filter <- result.table.list["k4"][[1]]
result.filter$halo.neuron <- result.filter$halo.astro <- 
  result.filter$halo.micro <- result.filter$halo.oligo <- 0
for(ii in seq(nrow(result.filter))){
  sample.id <- result.filter$sample_label[ii]
  sample.id.format <- paste0(unlist(strsplit(sample.id, "_"))[1:2], collapse = "_")
  which.halo.id <- halo.sample.prop$sample.id==sample.id.format
  halo.filter <- halo.sample.prop[which.halo.id,]
  result.filter$halo.neuron[ii] <- halo.filter$Excit+halo.filter$Inhib
  result.filter$halo.oligo[ii] <- halo.filter$Oligo
  result.filter$halo.astro[ii] <- halo.filter$Astro
  result.filter$halo.micro[ii] <- halo.filter$Micro
}
# get errors
result.filter$error.neuron <- result.filter$halo.neuron-result.filter$neuron
result.filter$abs.error.neuron <- abs(result.filter$error.neuron)
result.filter$error.oligo <- result.filter$halo.oligo-result.filter$oligo
result.filter$abs.error.oligo <- abs(result.filter$oligo)
result.filter$error.astro <- result.filter$halo.astro-result.filter$astro
result.filter$abs.error.astro <- abs(result.filter$error.astro)
result.filter$error.micro <- result.filter$halo.micro-result.filter$micro
result.filter$abs.error.micro <- abs(result.filter$error.micro)

summary(result.filter$abs.error.micro)

# append halo metadata
# total cells
# total cells per type

# append expression metadata
# total cells 
# total cells per type
