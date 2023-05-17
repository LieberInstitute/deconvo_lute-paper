#
# Get potential confounds from snRNAseq datasets
#
#

#----------------------
# prepare snrnaseq data
#----------------------
# load sce
sce <- get(load('./deconvo_method-paper/outputs/sce_prep_05-10-23_train.rda'))

#-------------------------------
# get counts and expressed genes
#-------------------------------
# total counts summaries

# expressed genes by sample

# save 

#------------------------------------
# get expression at top donor markers
#------------------------------------
# load donor marker data
list.markers.sn <- get(load("./deconvo_method-paper/outputs/list-markers-by-donor-k2-3-4_train.rda"))

# bulk marker expr
sample.id.index.vector <- seq(length(list.markers.sn))
list.sn.marker.expr <- lapply(sample.id.index.vector, function(index){
  sample.id <- names(list.markers.sn)[index]
  markers.vector <- unlist(list.markers.sn[[index]])
  sce.filter <- sce[rownames(sce) %in% markers.vector,sce$Sample==sample.id]
  lute:::.get_z_from_sce(sce.filter, "counts", "k2")
})
names(list.sn.marker.expr) <- paste0(sample.id.index.vector)
save(list.sn.marker.expr, file = "./deconvo_method-paper/outputs/list-donormarker-snexpr_dlpfc-train.rda")
