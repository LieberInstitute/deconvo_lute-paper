
#
# uses the following conditions
#
# * S : standard scale factors (neuron = 10, glial = 3)
# * scales, y and z : uses logcounts for y and z both
#
#

# dependencies
libv <- c("here", "lute", "dplyr", "ggplot2", "gridExtra", "SingleCellExperiment", "scuttle",
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
s.vector <- c("glial" = 12, "neuron" = 11,
              "glial_non_oligo" = 3, "astro" = 4, 
              "oligo" = 3, "micro" = 3)
sample.id.vector.sn <- unique(sce$Sample)
result.list <- lapply(sample.id.vector.sn, function(sample.id){
  sce.sample <- sce[,sce$Sample == sample.id]
  rse.sample <- rse[,rse$Sample == sample.id]
  # append logcounts
  sce.sample <- scuttle::logNormCounts(sce.sample)
  rse.sample <- scuttle::logNormCounts(rse.sample)
  assays(sce.sample)[["logcounts"]] <- as.matrix(assays(sce.sample)[["logcounts"]])
  assays(rse.sample)[["logcounts"]] <- as.matrix(assays(rse.sample)[["logcounts"]])
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
      y.expression <- assays(rse.markers)[["logcounts"]]
      colnames(y.expression) <- paste0(rse.markers$Sample, "_", 
                                       rse.markers$library_type, "_", 
                                       rse.markers$library_prep)
      
      # without s
      decon.no.scale <- lute(sce = sce.markers, y = y.expression,
                             celltype.variable = marker.type,
                             typemarker.algorithm = NULL,
                             assay.name = "logcounts")$deconvolution.results %>%
        as.data.frame()
      # with s
      decon.with.scale <- lute(sce = sce.markers, y = y.expression,
                               celltype.variable = marker.type, s = s.vector,
                               typemarker.algorithm = NULL,
                               assay.name = "logcounts")$deconvolution.results %>%
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

# result.table.name <- "deconvo-results-table_bulk-matched-sn-bydonor_cellsize-fract-1-0-9_train.rda"
result.table.name <- "deconvo-results-table_bulk-matched-sn-bydonor_s-10n-3g_scale-yz-logcounts_train.rda"
result.table.path <- paste0("./deconvo_method-paper/outputs/", result.table.name)
save(deconvo.results.list.k234, file = result.table.path)

#---------------
# map image data
#---------------
# load data
# deconvo results table
result.table.list <- get(load("./deconvo_method-paper/outputs/20_snrnaseq-bulk-matched_training/deconvo-results-table_bulk-matched-sn-bydonor_s-10n-3g_scale-yz-logcounts_train.rda"))
# image reference
halo.path <- "Human_DLPFC_Deconvolution/processed-data/03_HALO/halo_all.Rdata"
halo.all <- get(load(halo.path))
# filter halo.all
cell.types.filter <- c("Other")
halo.all <- halo.all[!halo.all$cell_type %in% cell.types.filter,]

# cell type proportions by sample
halo.sample.prop <- do.call(rbind, sapply(unique(halo.all$Sample), function(sample.id){
  table(halo.all[halo.all$Sample==sample.id,]$cell_type) %>% prop.table()
})) %>% as.data.frame()
halo.sample.prop$sample.id <- unique(halo.all$Sample)

# total cells by sample
neuron.types.vector <- c("Excit", "Inhib")
glial.types.vector <- c("Oligo", "OPC", "Micro", "Astro")
halo.cells.all <- table(halo.all$Sample) %>% as.data.frame()
halo.cells.neuron <- table(halo.all[halo.all$cell_type %in% neuron.types.vector,]$Sample) %>% as.data.frame()
halo.cells.glial <- table(halo.all[halo.all$cell_type %in% glial.types.vector,]$Sample) %>% as.data.frame()

# cell sizes by sample, type
halo.k2 <- halo.all
halo.k2$k2 <- ifelse(halo.k2$cell_type %in% neuron.types.vector, "neuron", "glial")
# aggregate cell size medians
halo.cellsize <- aggregate(halo.k2$Cell_Area, 
                           by = list(halo.k2$k2, halo.k2$Sample), FUN = median)
colnames(halo.cellsize) <- c("k2", "sample", "median.cellsize")
# jitterbox, by type
ggplot(halo.cellsize, aes(x = k2, y = median.cellsize)) + 
  geom_jitter(alpha = 0.4) + 
  geom_boxplot(alpha = 0, color = "cyan")
# jitterbox, by sample, type
ggplot(halo.cellsize, aes(x = k2, y = median.cellsize)) + 
  geom_jitter(alpha = 0.4) + 
  geom_boxplot(alpha = 0, color = "cyan") + facet_wrap(~sample)

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
  result.filter$halo.total.cells[ii] <- 
    halo.cells.all[halo.cells.all[,1]==sample.id.format,2]
  result.filter$halo.neuron.cells[ii] <- 
    halo.cells.neuron[halo.cells.neuron[,1]==sample.id.format,2]
  result.filter$halo.glial.cells[ii] <- 
    halo.cells.glial[halo.cells.glial[,1]==sample.id.format,2]
  result.filter$halo.glial.median.cellsize[ii] <- 
    halo.cellsize[halo.cellsize[,2]==sample.id.format & 
                    halo.cellsize[,1]=="glial", 3]
  result.filter$halo.neuron.median.cellsize[ii] <- 
    halo.cellsize[halo.cellsize[,2]==sample.id.format & 
                    halo.cellsize[,1]=="neuron", 3]
}
# get errors
result.filter$error.neuron <- result.filter$halo.neuron-result.filter$neuron
result.filter$abs.error.neuron <- abs(result.filter$error.neuron)
result.filter$error.glial <- result.filter$halo.glial-result.filter$glial
result.filter$abs.error.glial <- abs(result.filter$glial)
result.filter$cell.size.fraction <- 
  result.filter$halo.neuron.median.cellsize/
  result.filter$halo.glial.median.cellsize
# error summaries
median(result.filter[result.filter$scale==T,]$abs.error.neuron) # 0.3997985
median(result.filter[result.filter$scale==F,]$abs.error.neuron) # 0.1722251
median(result.filter[result.filter$scale==T,]$abs.error.glial) # 0.9148593
median(result.filter[result.filter$scale==F,]$abs.error.glial) # 0.3834312


# plots
# error by scale
ggplot(result.filter, aes(x = scale, y = abs.error.neuron)) + 
  geom_jitter(alpha = 0.4) + geom_boxplot(alpha = 0, color = "cyan")
# proportions, facet scale
ggplot(result.filter, aes(x = halo.neuron, y = neuron)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~scale)
# error by amount
ggplot(result.filter, aes(x = halo.neuron.cells, y = abs.error.neuron)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~scale) + geom_smooth(method = "glm")

# proportions, facet scale, color by condition
plot.data <- result.filter
sample.label.vector <- plot.data$sample_label
condition.vector <- sapply(sample.label.vector, function(sample.id){
  paste0(unlist(strsplit(sample.id, "_"))[3:4], collapse = "_")
})
plot.data$condition <- gsub("1", "", condition.vector)
ggplot(plot.data, aes(x = halo.neuron, y = neuron, color = condition)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~scale)

# error by cellsize fract
ggplot(result.filter, aes(x = cell.size.fraction, y = abs.error.neuron)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~scale) + geom_smooth(method = "glm")

# get scale error change by bulk sample
sample.id.vector <- unique(result.filter$sample_label)
sample.id.vector <- sample.id.vector[!grepl(".*Cyto1|.*Nuc1|.*Bulk1", sample.id.vector)]
plot.data <- do.call(rbind, lapply(sample.id.vector, function(sample.id){
  dff <- result.filter[grepl(sample.id, result.filter$sample_label),]
  value <- dff[dff$scale==FALSE,"abs.error.neuron"]-dff[dff$scale==TRUE,"abs.error.neuron"]
  c(sample.id, value, dff$halo.neuron.cells[1])
})) %>% as.data.frame()
colnames(plot.data) <- c("sample.id", "scale.error.diff", "halo.cells.neuron")
plot.data[,2] <- as.numeric(plot.data[,2])
plot.data[,3] <- as.numeric(plot.data[,3])
plot.data[,1] <- factor(plot.data[,1], levels = plot.data[,1][order(plot.data[,2])])

# barplot error diff
ggplot(plot.data, aes(x = sample.id, y = scale.error.diff)) + 
  geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# scatterplot error diff vs. num. halo neurons
ggplot(plot.data, aes(x = halo.cells.neuron, y = scale.error.diff)) +
  geom_hline(yintercept = 0) + geom_point(alpha = 0.4) + 
  geom_smooth(method = "glm")

ggplot(plot.data, aes(x = halo.cells.neuron, y = scale.error.diff)) +
  geom_hline(yintercept = 0) + geom_point(alpha = 0.4) + 
  geom_smooth(method = "glm") + geom_label(aes(label = sample.id))

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
