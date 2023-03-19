#!/usr/bin/env R

# Author: Sean Maden
#
# Get differentially expressed genes (DEGs) among bulk sample experiment conditions.
#

libv <- c("SingleCellExperiment", "SummarizedExperiment")
sapply(libv, library, character.only = T)

#----------
# load data
#----------
# get save directory path
save.path <- file.path("deconvo_method-paper", "outputs", 
                       "11_bulk-expression-analysis")

# load marker bulk expr
rsef.filename <- "rsef_k2-marker-expr_ro1-dlpfc.rda"
rsef <- get(load(file.path(save.path, rsef.filename)))

# load halo outputs
filename <- "halo_all.Rdata"
file.path <- file.path("Human_DLPFC_Deconvolution", "processed-data", "03_HALO", filename)
dfh <- get(load(path))

#-----------------------------------------------
# correlation of bulk and halo marker expression
#-----------------------------------------------
# prepare correlation matrix
dfcor <- data.frame(akt3_expr = dfp.med[dfp.med[,1]=="akt3_expr",4],
                    nuc_area = dfp.med[dfp.med[,1]=="nuc_area",4],
                    nuc_perim = dfp.med[dfp.med[,1]=="nuc_perim",4],
                    slide = dfp.med[dfp.med[,1]=="nuc_perim",3],
                    type = dfp.med[dfp.med[,1]=="nuc_perim",2])
# get plot object
lcor <- lapply(unique(dfcor$type), function(ti){
  mcor <- cor(dfcor[dfcor$type==ti,seq(3)], method = "spearman")
})

# make correlation heatmap
names(lcor) <- unique(dfcor$type)
lgg <- lapply(seq(3), function(ii){
  ggcorrplot(lcor[[ii]], type = "lower", lab = T, 
             title = names(lcor)[ii])
})
