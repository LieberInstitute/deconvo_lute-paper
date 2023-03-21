source("deconvo_method-paper/code/_halo-outputs-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
halo.output.table <- get(load(halo.output.path))
halo.output.table <- halo.output.table %>% as.data.frame()

# get transformed variables
# quantile transform marker counts by subject
samples.vector <- halo.output.table$Sample
unique.samples <- samples.vector %>% unique()
transformed.counts <- unlist(lapply(unique.samples, function(sample){
  message("working on sample: ", sample)
  filter <- halo.output.table$Sample==sample
  gene.marker.vector <- halo.output.table[filter, gene.marker.label]
  quantile.vector <-  gene.marker.vector %>% quantile(seq(0, 1, 1e-3))
  quantile.vector.names <- quantile.vector %>% names()
  quantile.vector.names <- gsub("%", "", quantile.vector.names)
  quantile.labels.vector <- quantile.vector.names %>% as.numeric()
  index.vector <- 1:(length(quantile.labels.vector)-1)
  quantile.labels.vector2 <- quantile.labels.vector[index.vector]
  new.vector <- rep("NA", length(gene.marker.vector))
  for(index in seq(length(quantile.labels.vector2))){
    which.marker <- gene.marker.vector >= quantile.labels.vector[index] & 
      gene.marker.vector < quantile.labels.vector[index+1]
    new.vector[which.marker] <- quantile.labels.vector[index+1]
  }
  return(new.vector)
}))
transformed.counts <- as.numeric(transformed.counts)
transformed.counts[is.na(transformed.counts)] <- median(transformed.counts, na.rm = T) # replace NAs
halo.output.table[,marker.quantile.variable] <- transformed.counts %>% as.numeric()
# get first transformation
halo.output.table[,transformed.marker.variable] <- transformation1(halo.output.table[,gene.marker.label])

# resave outputs table
save(halo.output.table, file = output.updated.path)





