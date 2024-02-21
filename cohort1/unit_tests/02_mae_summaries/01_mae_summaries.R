#!/usr/bin/env R
# Author: Sean Maden

#----------------
#
# load
#
#----------------

# all samples
mae.filename <- "mae_allsamples.rda"
mae.path <- file.path("outputs", "01_mae", mae.filename)
mae.all <- get(load(mae.path))

# removed samples
mae.filename <- "mae_analysis_append.rda"
mae.path <- file.path("outputs", "01_mae", mae.filename)
mae.analysis <- get(load(mae.path))

# training and validation sample ids
list.filename <- "list_snrnaseq_sampleid.rda"
list.file.path <- file.path("outputs", "00_preprocess", list.filename)
list.crossvalidate <- get(load(list.file.path))






#--------------------------------------------------
#
# test1, 
# snRNAseq samples count, no RNAscope filter
#
#--------------------------------------------------


#
assay.name <- "snrnaseq.k2.all"

all.samples.assay <- mae.all[[assay.name]]

analysis.samples.assay <- mae.analysis[[assay.name]]

snrnaseq.samples.vector <- unique(analysis.samples.assay$Sample)

training.samples.vector <- list.crossvalidate$train

validation.samples.vector <- list.crossvalidate$validation

# tests

# snrnaseq samples equals 17
length(snrnaseq.samples.vector)==17

# training samples from snrnaseq equals 12
length(training.samples.vector)==12

# validation samples from snrnaseq equals 12
length(validation.samples.vector)==5







#------------------------------------------
#
# test2, 
# snRNAseq and RNAscope sample count
#
#------------------------------------------

# assay name 1
assay.name1 <- "snrnaseq.k2.all"
all.samples.assay1 <- mae.all[[assay.name1]]
analysis.samples.assay1 <- mae.analysis[[assay.name1]]

# assay name 2
assay.name2 <- "sce.img"
all.samples.assay2 <- mae.all[[assay.name2]]
analysis.samples.assay2 <- mae.analysis[[assay.name2]]

###

length(unique(analysis.samples.assay2$Sample))==16

length(
  intersect(
    unique(analysis.samples.assay1$Sample), 
    unique(analysis.samples.assay2$Sample)
    )
  )

# tests

# training snrnaseq samples overlapping rnascope filtered snrnaseq samples equals 11
length(intersect(list.crossvalidate$train, all.samples.assay2$Sample))==11

# validation snrnaseq samples overlapping rnascope filtered snrnaseq samples equals 11
length(intersect(list.crossvalidate$validation, all.samples.assay2$Sample))==4






#------------------------------------------
#
# test3, 
# RNAscope count from cell.sizes
#
#------------------------------------------
assay.name1 <- "cell.sizes"

assay.name2 <- "snrnaseq.k2.all"

assay1 <- mae.analysis[[assay.name1]]

assay2 <- mae.analysis[[assay.name2]]

###

length(unique(assay1["sample_id",]))==16

length(
  intersect(
    
    unique(assay1["sample_id",]),
    
    unique(assay2[["Sample"]])
    
  )
)

#------------------------------------------
#
# test4, 
# Bulk samples condition summaries
#
#------------------------------------------

assay.name1 <- "bulk.rnaseq"

assay1 <- mae.analysis[[assay.name1]]

# total unique subjects
all.samples.vector.bulk.rnaseq <- unique(assay1$batch.id2)
length(all.samples.vector.bulk.rnaseq)==19

# total subjects per library prep
subjects.by.condition <- as.data.frame(
  table(assay1$library_prep, assay1$batch.id2))
sum(subjects.by.condition[subjects.by.condition[,1]=="Bulk",3])==38
sum(subjects.by.condition[subjects.by.condition[,1]=="Nuc",3])==37
sum(subjects.by.condition[subjects.by.condition[,1]=="Cyto",3])==38
# total subjects per library prep
subjects.by.condition <- as.data.frame(
  table(assay1$library_type, assay1$batch.id2))
sum(subjects.by.condition[subjects.by.condition[,1]=="polyA",3])==56
sum(subjects.by.condition[subjects.by.condition[,1]=="RiboZeroGold",3])==57

# sample source id summaries
# total sample source ids
length(unique(assay1$batch.id2))==19
# total sample source ids training
length(unique(assay1[,assay1$batch.id2 %in% training.samples.vector]$batch.id2))==11
# total sample source ids validation
length(
  unique(
    assay1[,assay1$batch.id2 %in% list.crossvalidate$validation]$batch.id2))==5


# by compartment
assay1$crossvalidation <- ifelse(assay1$batch.id2 %in% 
                                   training.samples.vector, "train", 
                                 ifelse(assay1$batch.id2 %in% 
                                          list.crossvalidate$validation, 
                                        "validate", "NA"))
# total
dft <- as.data.frame(table(assay1$library_prep, assay1$library_type))
for(index in seq(nrow(dft))){dft[index,]==19}
# crossvalidation sets

# train
filter1 <- 
  assay1$library_prep=="Bulk" & 
  assay1$library_type=="polyA" & 
  assay1$crossvalidation=="train"
length(
  unique(
    assay1[,filter1]$batch.id2))==11
dft <- as.data.frame(table(
  assay1$library_prep, 
  assay1$library_type, 
  assay1$crossvalidation))
dft[dft[,1]=="Bulk" & dft[,2]=="polyA" & dft[,3]=="train",4]==11
dft[dft[,1]=="Cyto" & dft[,2]=="polyA" & dft[,3]=="train",4]==11
dft[dft[,1]=="Nuc" & dft[,2]=="polyA" & dft[,3]=="train",4]==10
# validate
dft[dft[,1]=="Bulk" & dft[,2]=="polyA" & dft[,3]=="validate",4]==5
dft[dft[,1]=="Cyto" & dft[,2]=="polyA" & dft[,3]=="validate",4]==5
dft[dft[,1]=="Nuc" & dft[,2]=="polyA" & dft[,3]=="validate",4]==5
