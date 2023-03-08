#!/usr/bin/env R

# Author: Sean Maden
#
#

libv <- c("SummarizedExperiment")
sapply(libv, library, character.only = T)

#----------
# load data
#----------
# load bulk data
rse.filename <- "rse_gene.Rdata"
load.path <- file.path("Human_DLPFC_Deconvolution",
                       "processed-data",
                       "01_SPEAQeasy",
                       "round2_v40_2022-07-06",
                       "rse")
rse <- get(load(file.path(path, rse.filename)))

# get save directory path
save.path <- file.path("deconvo_method-paper", "outputs", 
                       "11_bulk-expression-analysis")

#-----------------------
# set analysis variables
#-----------------------
cd <- colData(rse)

table(cd$library_prep)
# Bulk Cyto  Nuc 
# 38   38   37 
table(cd$library_type)
# polyA RiboZeroGold 
# 56           57 
table(cd$library_type, cd$library_prep)
#               Bulk Cyto Nuc
# polyA          19   19  18
# RiboZeroGold   19   19  19

sample.variable <- paste0(cd$BrNum, "_", cd$location)
table(sample.variable)
# sample.variable
# Br2720_Mid Br2720_Post  Br2743_Ant  Br3942_Ant  Br3942_Mid  Br6423_Ant Br6423_Post  Br6432_Ant 
# 6           6           6           6           6           6           6           6 
# Br6432_Mid  Br6471_Ant  Br6471_Mid  Br6522_Mid Br6522_Post  Br8325_Ant  Br8325_Mid  Br8492_Mid 
# 6           6           6           6           6           6           5           6 
# Br8492_Post  Br8667_Ant  Br8667_Mid 
# 6           6           6
