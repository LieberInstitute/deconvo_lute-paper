library(SingleCellExperiment)

path <- "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030"

# source current functions
source.dpath <- "deconvo_method-paper/source"
for(fn in list.files(source.dpath)){source(file.path(source.dpath, fn))}

# load sce
sce.fpath <- "DLPFC_snRNAseq/processed-data/sce/sce_DLPFC.Rdata"
sce <- get(load(sce.fpath))

markerdata <- zsource_markerdata(sce)
typev <- zsource_type(zsource = sce, type.varname = "cellType_broad_hc")


markerdata <- DeconvoBuddies::get_mean_ratio2(sce, 
                                              cellType_col = "cellType_broad_hc", 
                                              assay_name = "logcounts",
                                              add_symbol = TRUE)
normcounts(sce) <- t(t(counts(sce))/sf)
logcounts(sce) <- log2(normcounts(sce)+1)
sce_assay <- as.matrix(SummarizedExperiment::assays(sce)[["logcounts"]])

# get new z
z.expt <- get_z_experiment(sce, method.markers = "mean_ratio", 
                           mr.assay = "logcounts", ngenes.byk = 25, 
                           type.varname = "cellType_broad_hc", 
                           summary.varname = "donor", k.summary.method = "mean",
                           z.summary.method = "mean", return.all = FALSE)