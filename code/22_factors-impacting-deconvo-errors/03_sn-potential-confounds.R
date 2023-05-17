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
library(scuttle)

# total counts summaries
# mean counts
count.summary <- summarizeAssayByGroup(sce, sce$Sample, statistics = "mean")
neuron.filter <- sce[["k2"]]=="neuron"
sce.filter <- sce[,neuron.filter]
neuron.count.summary <- summarizeAssayByGroup(sce.filter, sce.filter$Sample, statistics = "mean")
glial.filter <- sce[["k2"]]=="glial"
sce.filter <- sce[,glial.filter]
glial.count.summary <- summarizeAssayByGroup(sce.filter, sce.filter$Sample, statistics = "mean")
# total counts
df.counts <- data.frame(sample.id = colnames(assays(count.summary)[["mean"]]),
           total.counts = colSums(assays(count.summary)[["mean"]]),
           neuron.counts = colSums(assays(neuron.count.summary)[["mean"]]),
           glial.counts = colSums(assays(neuron.count.summary)[["mean"]]))

# expressed genes by sample
library(dplyr)
min.expr <- 5
df.expr.genes <- data.frame(sample.id = 
                              colnames(assays(count.summary)[["mean"]]),
                            total.expr.genes = 
                              apply(assays(count.summary)[["mean"]], 2, 
                                    function(ci){length(ci[ci>=min.expr])}) %>% as.numeric(),
                            neuron.expr.genes =
                              apply(assays(neuron.count.summary)[["mean"]], 2, 
                                    function(ci){length(ci[ci>=min.expr])}) %>% as.numeric(),
                            glial.expr.genes =
                              apply(assays(glial.count.summary)[["mean"]], 2, 
                                    function(ci){length(ci[ci>=min.expr])}) %>% as.numeric())

# save 
df.sn.confound <- cbind(df.expr.genes, df.counts) %>% as.data.frame()
save(df.sn.confound, file = 
       "./deconvo_method-paper/outputs/df-confounds-sn_dlpfc-train.rda")

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
