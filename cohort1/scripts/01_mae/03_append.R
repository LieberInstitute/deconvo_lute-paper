
#
# Append data to snRNAseq SingleCellExperiment objects
#


libv <- c("dplyr", "ggplot2")
sapply(libv, library, character.only = TRUE)

#-----
# load
#-----
mae.path <- "./outputs/01_mae/mae_analysis.rda"
mae <- get(load(mae.path))

# load snrnaseq filtered nucleus proportions
sn.path <- "./data/snRNA_cell_type_proportions.csv"
sn.proportions <- read.csv(sn.path)

#---------------------------
# true cell proportions list
#---------------------------

sample.id.vector <- unique(sn.proportions$Sample)

# k2 
list.df.true.k2 <- lapply(sample.id.vector, function(sample.id){
  sn.filter <- sn.proportions[sn.proportions$Sample==sample.id,]
  # sums
  total.neuron <- sum(sn.filter[sn.filter$cell_type %in% c("Excit", "Inhib"),]$n_cell_sn)
  total.glial <- sum(sn.filter[sn.filter$cell_type %in% c("Oligo", "Micro", "Astro"),]$n_cell_sn)
  total.cells <- sum(c(total.neuron, total.glial))
  # props
  prop.neuron <- total.neuron/total.cells
  prop.glial <- total.glial/total.cells
  # return
  dfr <- data.frame(glial = prop.glial, neuron = prop.neuron)
  rownames(dfr) <- "true_proportion"
  return(dfr)
})
names(list.df.true.k2) <- sample.id.vector
metadata(mae[["snrnaseq.k2.all"]])[["list.df.true.k2"]] <- list.df.true.k2

# k3
list.df.true.k3 <- lapply(sample.id.vector, function(sample.id){
  sn.filter <- sn.proportions[sn.proportions$Sample==sample.id,]
  # sums
  total.excit <- sum(sn.filter[sn.filter$cell_type %in% c("Excit"),]$n_cell_sn)
  total.inhib <- sum(sn.filter[sn.filter$cell_type %in% c("Inhib"),]$n_cell_sn)
  total.glial <- sum(sn.filter[sn.filter$cell_type %in% c("Oligo", "Micro", "Astro"),]$n_cell_sn)
  total.cells <- sum(c(total.excit, total.inhib, total.glial))
  # props
  prop.excit <- total.excit/total.cells
  prop.inhib <- total.inhib/total.cells
  prop.glial <- total.glial/total.cells
  # return
  dfr <- data.frame(glial = prop.glial, excit = prop.excit, inhib = prop.inhib)
  rownames(dfr) <- "true_proportion"
  return(dfr)
})
names(list.df.true.k3) <- sample.id.vector
metadata(mae[["snrnaseq.k3.all"]])[["list.df.true.k3"]] <- list.df.true.k3

# k4
list.df.true.k4 <- lapply(sample.id.vector, function(sample.id){
  sn.filter <- sn.proportions[sn.proportions$Sample==sample.id,]
  # sums
  total.excit <- sum(sn.filter[sn.filter$cell_type %in% c("Excit"),]$n_cell_sn)
  total.inhib <- sum(sn.filter[sn.filter$cell_type %in% c("Inhib"),]$n_cell_sn)
  total.oligo <- sum(sn.filter[sn.filter$cell_type %in% c("Oligo"),]$n_cell_sn)
  total.non.oligo.glial <- sum(sn.filter[sn.filter$cell_type %in% c("Micro", "Astro"),]$n_cell_sn)
  total.cells <- sum(c(total.excit, total.inhib, total.oligo, total.non.oligo.glial))
  # props
  prop.excit <- total.excit/total.cells
  prop.inhib <- total.inhib/total.cells
  prop.oligo <- total.oligo/total.cells
  prop.non.oligo.glial <- total.non.oligo.glial/total.cells
  # return
  dfr <- data.frame(oligo = prop.oligo, excit = prop.excit, inhib = prop.inhib, 
                    non.oligo.glial = prop.non.oligo.glial)
  rownames(dfr) <- "true_proportion"
  return(dfr)
})
names(list.df.true.k4) <- sample.id.vector
metadata(mae[["snrnaseq.k4.all"]])[["list.df.true.k4"]] <- list.df.true.k4

#------------------
# append cell sizes
#------------------
for(index in seq(3)){
  k.variable <- paste0("k", index+1)
  sce <- mae[[paste0("snrnaseq.",k.variable,".all")]]
  sample.id.vector <- unique(sce[["Sample"]])
  df.sizes.all <- do.call(rbind, lapply(sample.id.vector, function(sample.id){
    scef <- sce[,sce[["Sample"]]==sample.id]
    cell.types <- unique(scef[[k.variable]])
    df.sizes <- do.call(rbind, lapply(cell.types, function(cell.type){
      counts <- assays(
        scef[,scef[[k.variable]]==cell.type])[["counts"]]
      mean(colSums(counts))
    }))
    df.sizes <- as.data.frame(df.sizes)
    colnames(df.sizes) <- "size"
    df.sizes$k4 <- cell.types
    df.sizes$sample.id <- sample.id
    df.sizes$size.type <- "mean_library_counts"
    return(df.sizes)
  }))
  df.sizes.all <- as.data.frame(df.sizes.all)
  df.sizes.all$ktype <- k.variable
  assay.name <- paste0("snrnaseq.",k.variable,".all")
  metadata(mae[[assay.name]])[["cell.sizes"]] <- df.sizes.all
}

#-----
# save
#-----
mae.new.path <- "./outputs/01_mae/mae_analysis_append.rda"
save(mae, file = mae.new.path)
