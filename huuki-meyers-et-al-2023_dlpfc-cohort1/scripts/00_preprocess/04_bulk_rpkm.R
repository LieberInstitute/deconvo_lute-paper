# make rpkm object
# load rpkm counts
rpkm.filename <- "rpkmCounts_Human_DLPFC_Deconvolution_n113.rda"
rpkm.path <- file.path("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/Human_DLPFC_Deconvolution/processed-data/01_SPEAQeasy/round2_v25_2023-04-05/count_objects", rpkm.filename)
file.exists(rpkm.path)
load(rpkm.path)
# make rse from rpkm gene counts
libv <- c("SummarizedExperiment", "dplyr"); sapply(libv, library, character.only = T)
rse.rpkm <- SummarizedExperiment(assays = geneRpkm, rowData = geneMap) %>% as("RangedSummarizedExperiment")
# save rse from rpkm gene counts

rse.rpkm.filename.save <- "rse-rpkmCounts_Human_DLPFC_Deconvolution_n113.rda"
rse.rpkm.save.newpath <- file.path("users/smaden/", rse.rpkm.filename.save)
save(rse.rpkm, file = rse.rpkm.save.newpath)
