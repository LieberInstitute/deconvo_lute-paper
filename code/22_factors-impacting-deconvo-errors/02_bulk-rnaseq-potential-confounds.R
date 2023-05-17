
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

#--------------------------------------------
# samples with fewest counts, expressed genes
#--------------------------------------------