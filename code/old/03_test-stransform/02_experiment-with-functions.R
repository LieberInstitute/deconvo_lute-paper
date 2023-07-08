#!/usr/bin/env R

#
# Running the stransformation experiment using more functions, updated pseudobulk data, etc.
# 

#------
# paths
#------
proj.dpath <- "deconvo_method-paper"
save.dpath <- file.path(proj.dpath, "outputs/03_test-stransform")
source.dpath <- file.path(proj.dpath, "source")
script.fnamev <- c("z_methods.R", "z_transform.R", 
                   "make_example_data.R", 
                   "z_figures.R", "y_methods.R")

#-------
# params
#-------
celltype.treg.varname <- "celltype.treg"

#-----
# load
#-----
# source key functions
for(scripti in script.fnamev){source(file.path(source.dpath, scripti))}

#------------------
# make pb expt data
#------------------
# set number of marker genes
nfeatures <- 1000
# get counts matrix
scale.cell <- 1000; nk <- 4
num.cells.total <- nk*scale.cell
countsv <- sample(seq(100), nfeatures*num.cells.total, replace = T)
ct <- matrix(countsv, nrow = nfeatures) # counts matrix
# make sef
sef <- SummarizedExperiment(assays = list(counts = ct))
# define the cell labels
celltype.varname <- "celltype"
sef[[celltype.varname]] <- rep(c("excit", "inhib", "oligo", "other"),
                         each = scale.cell)

#----------------------
# get pseudobulk series
#----------------------
# make sample series
datv.j1 <- c(10,10,2,5) # skew towards neurons
datv.j2 <- c(20, 10, 2, 5) # more skew towards excit than inhib
datv.j3 <- c(5, 5, 5, 5) # all equal
lpb <- get_lpb(sef, datv = c(datv.j1, datv.j2, datv.j3), 
               ctvarname = celltype.varname)

#--------------
# get pb report
#--------------
z <- lpb$pi_pb











