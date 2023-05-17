
#
# Checks potential confounding variables from bulk RNA-seq data
# * total counts
# * total expressed genes
#

#-----------------------------
# prepare bulk rnaseq rse data
#-----------------------------
# bulk rnaseq
# paths
rse.path <- "Human_DLPFC_Deconvolution/processed-data/rse/rse_gene.Rdata"
rse <- get(load(rse.path))

#-------------------------------
# get counts and expressed genes
#-------------------------------
# total counts by sample
total.counts <- data.frame(sample.id = colnames(rse),
                           total.counts = colSums(assays(rse)[["counts"]]))

# expressed genes by sample
get_expressed_genes <- function(rse, min.expression = 5){
  m <- apply(assays(rse)[["counts"]], 2, function(ci){
    length(ci[ci>=min.expression])
  })
  df <- data.frame(sample.id = names(m), expressed.genes = as.numeric(m))
}
expressed.genes <- get_expressed_genes(rse)

# save 
list.bulk.confounds <- list(total.counts = total.counts,
                            expressed.genes = expressed.genes)
save(list.bulk.confounds, 
     file = "./deconvo_method-paper/outputs/list-bulk-confounds_dlpfc-train.rda")

#------------------------------------
# get expression at top donor markers
#------------------------------------
# load donor marker data
list.markers.sn <- get(load("./deconvo_method-paper/outputs/list-markers-by-donor-k2-3-4_train.rda"))

# bulk marker expr
sample.id.index.vector <- seq(length(list.markers.sn))
list.bulk.marker.expr <- lapply(sample.id.index.vector, function(index){
  sample.id <- names(list.markers.sn)[index]
  markers.vector <- unlist(list.markers.sn[[index]])
  rse[rownames(rse) %in% markers.vector,]
})
names(list.bulk.marker.expr) <- paste0(sample.id.index.vector)
save(list.bulk.marker.expr, file = "./deconvo_method-paper/outputs/list-donormarker-bulkexpr_dlpfc-train.rda")
