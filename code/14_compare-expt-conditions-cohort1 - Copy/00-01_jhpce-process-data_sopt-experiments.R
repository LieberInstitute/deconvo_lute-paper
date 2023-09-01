#!/usr/bin/env R

# Author: Sean Maden 
# 
# Run experiment script on JHPCE. Get de novo (i.e. not batch-robust) markers for given K.
# 
# On JHPCE, run:
# qrsh -l h_vmem=100G,mem_free=100G
# module load conda_R/4.3
# R
#
#

libv <- c("devtools", "scuttle")
sapply(libv, library, character.only = T)
devtools::install_github("metamaden/lute")
library(lute)

#---------------------------
# prep experiment parameters
#---------------------------
# params
rd.geneid.var.y <- "Symbol"
validation.sample.id <- c("Br6432", "Br6522", "Br8667")
assay.name.lute <- "logcounts"
group.id.variable <- group.name.experiment <- "sample.id"
validation.sample.id <- c("Br6432", "Br6522", "Br8667")
rd.geneid.var.y <- "Symbol"
rd.geneid.var.sce <- "gene_name"
y.group.variable.name <- "batch.id2"

#-----
# load
#-----
project.filepath <- "dcs04/lieber/lcolladotor/deconvolution_LIBD4030/"
user.filepath <- "users/smaden"

# get assay data
# load sce data
sce.path <- "dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/sce/sce_DLPFC.Rdata" # file.path(project.filepath, "")
sce <- get(load(sce.path))
# append sce k2
celltypevar <- "cellType_broad_hc"
sce[["k2"]] <- ifelse(grepl("^Excit.*|^Inhib.*", sce[[celltypevar]]), 
                      "neuron", "glial")
sce <- sce[,!sce[["k2"]]=="other"]
dim(sce)

# load bulk data
# y.unadj <- mae[["bulk.rnaseq"]]
y.path <- "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/Human_DLPFC_Deconvolution/processed-data/01_SPEAQeasy/round2_v25_2023-04-05/count_objects/rse_gene_Human_DLPFC_Deconvolution_n113.Rdata" # file.path(project.filepath, "")
y.unadj <- get(load(y.path))
rownames(y.unadj) <- rowData(y.unadj)[,rd.geneid.var.y]
y.unadj <- y.unadj[rownames(y.unadj) %in% rownames(sce),]
y.unadj <- y.unadj[order(match(rownames(y.unadj), rownames(sce))),]
y.unadj <- logNormCounts(y.unadj)

# format y sample.id
cd.y <- colData(y.unadj)
head(cd.y[,seq(5)])
y.id.vector <- cd.y[,1]
y.sample.id <- unlist(strsplit(y.id.vector, "_"))[2]
y.region <- tolower(unlist(strsplit(y.id.vector, "_"))[3])
colData(y.unadj)$sample.id <- paste0(y.sample.id, "_", y.region)
colData(y.unadj)$BrNum <- y.sample.id
# y.unadj$sample.id <- y.unadj$batch.id2
y.train <- y.unadj[,!y.unadj$BrNum %in% validation.sample.id]
y.validate <- y.unadj[,y.unadj$BrNum %in% validation.sample.id]
# load rnascope data
# df.rnascope.kdata
df.rn.path <- file.path(user.filepath, "df-rnascope-info_cohort1.rda")
df.rn <- get(load(df.rn.path))

# source scripts
script.names <- c("00_source_example-soptimize-framework-lute.R",
                  "00_source-deconvo-plots.R",
                  "01_source-example-crossval-sopt.R")
script.paths <- file.path(user.filepath, script.names)
sapply(script.paths, source)
# source functions
# source deconvo buddies info

colData(sce)[[group.id.variable]] <- sce[["BrNum"]]
colData(y.unadj)[[group.id.variable]] <- y.unadj[["BrNum"]]
# sample.id.vector <- unique(sce[[group.id.variable]])

#------------------
# 1. K2 Experiments
#------------------
# Prep k2 parameters
# get list.df.true, rnascope data
celltype.variable <- k.variable.name <- "k2"
unique.types <- unique(sce[[celltype.variable]])
list.df.true <- df.true.list(df.rn, unique(y.unadj$group.id), celltype.variable, unique.types)

# 1C. Experiment: 
# * P_true : RNAscope proportions
# * G : De novo K2 markers (N = 20 genes/celltype)
# * Z : snRNAseq matched or unmatched

k2.1c.result.list <- kmatch_experiment(k.variable.name = k.variable.name,
                                        sce = sce, 
                                        sample.id.vector = sample.id.vector, 
                                        list.df.true = list.df.true, 
                                        y.eset = y.unadj, y.train = y.train, 
                                        y.validate = y.validate, 
                                        assay.name = assay.name.lute, 
                                        group.name = group.name.experiment)

k2.1c.result.path <- file.path(user.path, "k2_1c_result_list.rda")
save(k2.1c.result.list, file = k2.1c.result.path)

# 1D. Experiment: 
# * P_true : snRNAseq proportions
# * G : De novo K2 markers (N = 20 genes/celltype)
# * Z : snRNAseq matched or unmatched

k2.1d.result.list <- kmatch_experiment(k.variable.name = k.variable.name,
                                       sce = sce, 
                                       sample.id.vector = sample.id.vector, 
                                       list.df.true = list.df.true, 
                                       y.eset = y.unadj, y.train = y.train, 
                                       y.validate = y.validate, 
                                       assay.name = assay.name.lute, 
                                       celltype.variable = k.variable.name, 
                                       group.name = group.name.experiment)

k2.1d.result.path <- file.path(user.path, "k2_1d_result_list.rda")
save(k2.1d.result.list, file = k2.1d.result.path)

# save results
k2.1c.res.path <- ""
k2.1d.res.path <- ""
save(k2.1c.result.list, file = k2.1c.res.path)
save(k2.1d.result.list, file = k2.1c.res.path)

#------------------
# 1. K3 Experiments
#------------------
