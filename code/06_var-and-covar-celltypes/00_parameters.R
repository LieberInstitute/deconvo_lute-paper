# dependencies
libv <- c("here", "lute", "dplyr", "ggplot2", "gridExtra", "SingleCellExperiment", 
          "SummarizedExperiment", "scran", "MultiAssayExperiment", "ComplexHeatmap", "pheatmap")
sapply(libv, library, character.only = TRUE)

# save path
prepped.data.path <- here("deconvo_method-paper", "outputs", "01_prepare-datasets")
save.path <- here("deconvo_method-paper", "outputs", "05_cell-count-versus-cell-size")

#-------------------
# dlpfc cohort1 data
#-------------------
# multi assay experiment path
mae.filename <- "mae_final.rda"
mae.path <- here("deconvo_method-paper", "outputs", "01_prepare-datasets", mae.filename)
# dlpfc markers path
sce.markers.list.path <- here("deconvo_method-paper", "outputs", "01_prepare-datasets", "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda")

#-------------------
# dlpfc cohort2 data
#-------------------
# mrb sce path
sce.mrb.path <- here("deconvo_method-paper", "outputs", "01_prepare-datasets", "sce-mrb_dlpfc.rda")

#-----------------
# helper functions
#-----------------
z_list <- function(sce, s.vector, assay.name = "logcounts", celltype.variable = "k2"){
  z1 <- lute:::.get_z_from_sce(sce, assay.name, celltype.variable)
  z2 <- lute:::.zstransform(z1, s.vector)
  list.return <- list("unadj" = z1, "adj" = z2)
  return(list.return)
}

z_list_bydonor <- function(sce, sample.id.variable = "Sample", 
                           s.vector = c("glial" = 3, "neuron" = 10), 
                           assay.name = "logcounts", celltype.variable = "k2"){
  unique.sample.id.vector <- unique(sce[[sample.id.variable]])
  zlist.bydonor <- lapply(unique.sample.id.vector, function(sample.id){
    z_list(sce[,sce[[sample.id.variable]]==sample.id], s.vector, assay.name, celltype.variable)
  })
  names(zlist.bydonor) <- unique.sample.id.vector
  return(zlist.bydonor)
}

covcor_from_zlist <- function(zlist, type = "cov"){
  if(is(zlist[[1]], "matrix")){
    lcov <- lapply(zlist, function(item){cov(item)})
    lcor <- lapply(zlist, function(item){cor(item)})
    lcov.diff.unadj.minus.adj <- lcov$unadj - lcov$adj
    lcor.diff.unadj.minus.adj <- lcor$unadj - lcor$adj
  } else{
    lcov <- lapply(zlist.bydonor, function(item){
      list("unadj" = cov(item[[1]]), 
           "adj" = cov(item[[2]]))})
    lcor <- lapply(zlist.bydonor, function(item){
      list("unadj" = cor(item[[1]]), "adj" = cor(item[[2]]))})
    lcov.diff.unadj.minus.adj <- lapply(zlist.bydonor, function(item){
      cov(item[[1]])-cov(item[[2]])})
    lcor.diff.unadj.minus.adj <- lapply(zlist.bydonor, function(item){
      cor(item[[1]])-cor(item[[2]])})
    names(lcov) <- names(lcor) <-
      names(lcov.diff.unadj.minus.adj) <-
      names(lcor.diff.unadj.minus.adj) <- names(zlist.bydonor)
  }
  list(lcov = lcov, lcor = lcor, 
       lcov.diff.unadj.minus.adj = lcov.diff.unadj.minus.adj,
       lcor.diff.unadj.minus.adj = lcor.diff.unadj.minus.adj)
}


