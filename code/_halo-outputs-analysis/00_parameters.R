#-------
# header
#------- 

libv <- c("nlme", "ggplot2", "gridExtra", "dplyr")
sapply(libv, library, character.only = TRUE)

# set the save directory
save.path <- here("deconvo_method-paper", "outputs", "09_manuscript")

# set the halo data path
halo.output.file.name <- "halo_all.Rdata"
halo.output.path <- here("Human_DLPFC_Deconvolution", "processed-data", "03_HALO", 
                         halo.output.file.name)

# cell labels
labels <- c("Endo" = "CLDN5", "Astro" = "GFAP", "Inhib" = "GAD1", "Excit" = "SLC17A7", 
            "Micro" = "TMEM119", "Oligo" = "OLIG2")

# helper functions
quantile_normalization1 <- function(levels.vector){
  # get quantiles by level, e.g. quantiles by sample (a.k.a. "slide")
  unique.levels <- levels.vector %>% unique()
  transformed.counts <- unlist(lapply(unique.levels, function(level){
    message("working on sample: ", level)
    filter <- levels.vector==level
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
  transformed.counts[is.na(transformed.counts)] <- median(transformed.counts, na.rm = T)
  transformed.counts <- transformed.counts %>% as.numeric()
  return(transformed.counts)
}

quantile_normalization2 <- function(table){
  require(preprocessCore)
  df_rank <- table %>% apply(2, rank, ties.method = "min")
  df_sorted <- table %>% apply(2, sort) %>% data.frame()
  df_mean <- df_sorted %>% apply(1, mean)
  index_to_mean <- function(my_index, my_mean){return(my_mean[my_index])}
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean = df_mean)
  
  
  rownames(df_final) <- rownames(df)
  return(df_final)
}

inverse_maximum_difference_transformation <- function(variable.vector){
  # marker transformation 1 : 1 / (max + 1) - value
  1/(max(variable.vector + 1) - variable.vector)
}

transformation2 <- function(variable){}

#--------
# scripts
#--------

# 01
sample.id.label <- "Sample"
cell.area.variable <- "Nucleus_Area"
gene.marker.label <- "AKT3_Copies"
#
output.updated.filename <- "halo_updated_path.Rdata"
output.updated.path <- here(save.path, output.updated.filename)
#
cell.area.log.variable <- "log10_nucleus_area"
marker.quantile.variable <- "akt3.copies.quantile.scale"
transformed.marker.variable <- "akt3.transformation1"
halo.quantiles.jpg.file.name <- ""

# 02

# 99 quantile scale summaries






