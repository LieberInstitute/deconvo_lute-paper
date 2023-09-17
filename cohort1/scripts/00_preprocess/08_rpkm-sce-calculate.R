library(edgeR)


# multi assay experiment path
mae.filename <- "mae_with-rpkm_additional-data_final.rda"
mae.path <- file.path("deconvo_method-paper", "outputs", "01_prepare-datasets", mae.filename)
mae <- get(load(mae.path))

# load data
sce <- mae[["sn1.rnaseq"]]
rse <- mae[["bulk.rnaseq"]]; rse <- rse[rownames(sce),]
sce.counts <- assays(sce)[["counts"]]
# load library sizes from complete snrnaseq object
sce.libsize.path <- "./deconvo_method-paper/outputs/01_prepare-datasets/total-counts-library-sizes_sce-dlpfc-cohort1.rda"
sce.libsize <- get(load(sce.libsize.path))
# set params
marker.lengths <- rowData(rse)$Length
library.sizes <- rep(1000, ncol(sce.counts))

#--------------
# full sce data
#--------------
# match libsize data
ls.final <- sce.libsize[colnames(sce.counts)]
identical(names(ls.final), colnames(sce.counts))
# get input
sce.input <- DGEList(counts = sce.counts, lib.size = ls.final,
                     genes = data.frame(Length=rowData(rse)$Length))
# get rpkm
sce.rpkm <- rpkm(y = sce.input, gene.length = sce.input$genes$Length, log = T)

# save
sce.name <- "sce-rpkm-counts_dlpfc-cohort1.rda"
sce.path <- file.path("./deconvo_method-paper/outputs/01_prepare-datasets", sce.name)
save(sce.rpkm, file = sce.path)

#-----------------------------------
# try with z reference dataset -- k2
#-----------------------------------
# library size by k2 type
libsize.vector <- sapply(unique(sce[["k2"]]), function(cell.type){
  celltype.filter <- sce[["k2"]]==cell.type
  mean(sce.libsize[celltype.filter])})
sce.z.counts <- lute:::.get_z_from_sce(sce, "counts", "k2")
rpkm.input <- DGEList(counts = sce.z.counts, lib.sizes = libsize.vector,
                      genes = data.frame(Length=marker.lengths))
sce.z.rpkm <- rpkm(y = rpkm.input, gene.length = rpkm.input$genes$Length, log = T)
# save
sce.z.rpkm.path <- file.path("./deconvo_method-paper/outputs/01_prepare-datasets/sce-z-rpkm-k2_dlpfc-cohort1.rda")
save(sce.z.rpkm, file = sce.z.rpkm.path)

# check bulk rpkm
# assays(mae[["bulk.rpkm.rnaseq"]])[[1]]

#-----------------------------------
# try with z reference dataset -- k3
#-----------------------------------

#-----------------------------------
# try with z reference dataset -- k4
#-----------------------------------
