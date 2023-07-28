library(scran)
sce.path <- "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/sce/sce_DLPFC.Rdata"
save.path <- "/users/smaden/sce-dlpfc-cohort1_qc_metrics.rda"
sce <- get(load(sce.path))
metrics <- perCellQCMetrics(sce)
save(metrics, file = save.path)