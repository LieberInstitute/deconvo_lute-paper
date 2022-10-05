#!/usr/bin/env R

#
# Setting up the S-transformation experiment. We wish to compare the pi_est in
# 3 conditions:
# 1. no S-transformation
# 2. "static" S-transformation (e.g. using cell type means)
# 3. "randomized" S-transformation (e.g. sampling rnorm by cell type)
# 
#

#------
# paths
#------
proj.dpath <- "deconvo_method-paper"
source.dpath <- file.path(proj.dpath, "source")
save.dpath <- file.path(proj.dpath, "outputs/03_test-stransform")

lz.fname <- "lz-mr_expt-stransform_dlpfc-ro1.rda"
script.fnamev <- c("z_methods.R", "z_transform")

sce.fpath <- "DLPFC_snRNAseq/processed-data/sce/sce_DLPFC.Rdata"

#-------
# params
#-------
celltype.varname <- "cellType_broad_hc"

#-----
# load
#-----
# source key functions
for(scripti in script.fnamev){source(file.path(source.dpath, scripti))}

# single cell experiment
sce <- get(load(sce.fpath))

#-------------------------------------------------
# make new cell type -- make types per tregs paper
#-------------------------------------------------
varv <- sce[[celltype.varname]]
table(varv)
#  Astro EndoMural     Micro     Oligo       OPC     Excit     Inhib
#   3979      2157      1601     32051      1940     24809     11067

new.ct <- rep("other", ncol(sce))
net.ct[which(sce[[celltype.varname]])]


#-----------------------------------------------------
# make lz -- list containing z table and other outputs
#-----------------------------------------------------

# get detailed zexpt data
lz <- get_z_experiment(sce, method.markers = "mean_ratio", 
                       mr.assay = "logcounts", ngenes.byk = 25, 
                       type.varname = "cellType_broad_hc", 
                       summary.varname = "BrNum", k.summary.method = "mean",
                       z.summary.method = "mean", return.all = TRUE,
                       save.dpath = save.dpath)



