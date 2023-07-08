libv <- c("here", "dplyr", "biomaRt")
sapply(libv, library, character.only = TRUE)

# sutton markers
load.name <- "sigsBrain.rda"
load.path <- here("deconvo_method-paper", "outputs", 
                  "15_k2-simulations_within-sample-matched", load.name)
sutton <- get(load(load.path))
# format sutton data
glial.vector <- c("Astrocytes", "Microglia", "Oligodendrocytes", "OPCs")
neuron.vector <- c("Excitatory", "Inhibitory")
mart<- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL"))
index.vector <- seq(length(sutton))
sutton.format <- lapply(index.vector, function(index){
  data <- sutton[[index]]
  name <- names(sutton)[index]
  message("working on dataset: ", name)
  # get k2 types
  cond.glial <- colnames(data) %in% glial.vector
  data$glial <- rowMeans(as.matrix(data[,cond.glial,drop=F])) + 1
  data$neuron <- as.matrix(data[,"Neurons"]) + 1
  cond.neuron <- intersect(colnames(data), neuron.vector) %>% length()
  if(cond.neuron==2){
    data$neuron <- rowMeans(as.matrix(data[,neuron.vector])) + 1
  }
  data$mean.ratio.over.neuron <- data$Neurons/data$glial
  data$mean.ratio.over.glial <- data$glial/data$Neuron
  # convert gene names
  ensemblGenes <- rownames(data)
  ens <- getBM(filters="ensembl_gene_id", 
               attributes=c("ensembl_gene_id", "hgnc_symbol"),
               values=ensemblGenes,
               mart=mart)
  ens <- ens[!duplicated(ens[,2]),]
  filter.data <- rownames(data) %in% ens$ensembl_gene_id
  data <- data[filter.data,]
  data <- data[order(match(rownames(data), ens$ensembl_gene_id)),]
  cond <- identical(rownames(data), ens$ensembl_gene_id)
  if(cond){rownames(data) <- ens$hgnc_symbol}
  return(data)
})
names(sutton.format) <- names(sutton)
# save
format.name <- "sutton-brain-markers-formatted.rda"
format.path <- here("deconvo_method-paper", "outputs", 
                    "18_updated-marker-sets-analyses", format.name)
save(sutton.format, file = format.path)
