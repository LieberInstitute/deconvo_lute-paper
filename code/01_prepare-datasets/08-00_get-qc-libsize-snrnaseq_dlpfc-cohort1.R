library(scran)
sce.path <- "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/sce/sce_DLPFC.Rdata"
save.path <- "/users/smaden/total-counts-library-sizes_sce-dlpfc-cohort1.rda"
sce <- get(load(sce.path))
# metrics <- perCellQCMetrics(sce)
metrics <- colSums(counts(sce))
names(metrics) <- colnames(sce)
save(metrics, file = save.path)