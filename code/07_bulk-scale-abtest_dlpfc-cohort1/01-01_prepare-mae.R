#!/usr/bin/env R

#
#
#


# load mae 
#new.mae.filename <- "mae_with-rpkm_additional-data_final.rda"
#mae.final.filepath <- here("deconvo_method-paper", "outputs", "01_prepare-datasets", new.mae.filename)
#mae <- get(load(mae.final.filepath))

# experiment variables
assay.name <- "logcounts"
deconvolution.algorithm <- "nnls"
cell.type.variable <- "k2"

#-----------
# unpack mae
#-----------
# get bulk expression
rse.counts <- mae[["bulk.rnaseq"]]
rse.rpkm <- mae[["bulk.rpkm.rnaseq"]]
names(assays(rse.counts)) <- names(assays(rse.rpkm)) <- "counts"

# snrnaseq reference -- using same reference across experiments
sce.iter <- mae[["sn1.rnaseq"]]
sce.iter <- logNormCounts(sce.iter)

# prep rnascope data
# get true proportions from rnascope data
rnascope <- mae[["rnascope.image"]]
# append k2 label
cd <- colData(rnascope)
cd$k2 <- "NA"
cd$k2 <- ifelse(grepl("Excit|Inhib", cd$cell_type), "neuron",
                ifelse(grepl("Endo|Oligo|Micro", cd$cell_type), "glial", "NA"))
colData(rnascope) <- cd
# filter na
rnascope <- rnascope[,!rnascope$k2=="NA"]
sample.id.vector <- unique(rnascope$Sample)
# make rnscope data (ROWS = CELLTYPE;SAMPLE_ID)
df.rn <- do.call(rbind, lapply(sample.id.vector, function(sample.id){
  rnf <- rnascope[,rnascope$Sample==sample.id]
  # proportions
  df.prop <- table(rnf[[cell.type.variable]], rnf$Sample) %>% prop.table()
  # sizes
  df.size <- aggregate(data.frame(area = assays(rnf)[["Nucleus_Area"]][1,]), 
                       by = list(cell_type = rnf[[cell.type.variable]]), FUN = "median")
  df.iter <- cbind(df.prop, df.size)
  df.iter <- df.iter[,c(1,2,3,5)]
  colnames(df.iter) <- c("cell_type", "sample_id", "true_proportion", "cell_size")
  df.iter
}))
