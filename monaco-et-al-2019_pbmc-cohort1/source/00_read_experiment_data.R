#!/usr/bin/env R

# Author: Sean Maden
#
# Read in the experiment elements.
#


processZ <- function(path){
  # processZ
  #
  # Get reference from file.
  
  # load
  zref <- 
    read.csv(
      "./data/monaco_et_al_2019/manuscript/abisseq_rnaseq_cell_types_references_k17.csv")
  
  # process
  zref <- zref[!zref[,1]=="12:00 AM",]
  rownames(zref) <- zref[,1]
  zref <- zref[,c(2:ncol(zref))]
  
  return(zref)
}

processP <- function(path, zFilter){
  # processP
  # Get proportions from file.
  #
  
  p <- read.csv(path)
  
  # process
  # set sample labels as row names
  rownames(p) <- p[,1]
  p <- p[,seq(2,ncol(p))]
  # subset cell types
  p <- p[,colnames(p) %in% colnames(zFilter)]
  # convert fractions
  sampleIdVector <- rownames(p)
  cellTypeVector <- colnames(p)
  p <- apply(p, 1, function(ri){ri/sum(ri)}) |> as.data.frame()
  
  return(p)
}

processY <- function(path, assayName = "tpm"){
  # processY
  #
  # Get bulk expression from file.
  
  # load
  yTable <- read.table(path)
  
  
  
  #---------
  
  
  # process
  
  
  #---------
  assayList <- list(yTable)
  names(assayList) <- assayName
  ySe <- SummarizedExperiment(assays = assayList)
  returnList <- list(
    yTable=yTable,
    ySe=ySe
  )
  
  # get coldata
  # get s13 phenotype info from column labels
  phenoDataS13 <- colnames(yTable)
  table(gsub(".*_", "", phenoDataS13))
  phenoDataS13 <- data.frame(sample.id = phenoDataS13,
                             source.id = gsub("_.*", "", phenoDataS13),
                             sample.type = gsub(".*_", "", phenoDataS13))
  phenoDataS13$tissue.type <- 
    ifelse(phenoDataS13$sample.type=="PBMC", "PBMC", "immune_cell")
  phenoDataS13$tissue.type.detail <- 
    ifelse(phenoDataS13$sample.type=="PBMC", "PBMC", phenoDataS13$sample.type)
  # append summary statistics to pheno data
  yTable <- as.matrix(yTable)
  phenoDataS13$library.size <- colSums(yTable)
  phenoDataS13$mean.expression <- colMeans(yTable)
  phenoDataS13$median.expression <- colMedians(yTable)
  phenoDataS13$sd.expression <- colSds(yTable)
  phenoDataS13$num.na.expression <- colAnyNAs(yTable)
  phenoDataS13$num.zero.expression <- unlist(apply(yTable, 2, function(ci){length(ci[ci==0])}))
  rownames(phenoDataS13) <- phenoDataS13[,1]
  # append pheno data to colData for SummarizedExperiment
  colData(ySe) <- DataFrame(phenoDataS13)
  
  
  return(returnList)
}

yMapMarkers <- function(bulkSummarizedExperiment){
  # yMapMarkers
  #
  # Get gene symbols for bulk SummarizedExperiment with ensembl gene ids.
  #
  #
  
  require(biomaRt)
  
  
  # begin rowdata
  assayData <- assays(bulkSummarizedExperiment)[[1]]
  newRowData <- data.frame(ensembl_gene_id_version = rownames(assayData))
  rownames(newRowData) <- newRowData[,1]
  rowData(bulkSummarizedExperiment) <- DataFrame(newRowData)
  
  # map gene ids to symbols
  ensemblMart <- useEnsembl(
    biomart = "genes", dataset = "hsapiens_gene_ensembl")
  geneIdVector <- rownames(bulkSummarizedExperiment)
  geneIdVector <- gsub("\\..*", "", geneIdVector)
  rowDataMaps <- getBM(
    attributes = c('ensembl_gene_id', 'hgnc_symbol'),
    filters = 'ensembl_gene_id',
    values = as.character(geneIdVector), 
    mart = ensemblMart
  )
  
  # map symbols for rowdata
  rowDataSymbol <- rowDataMaps[,2]
  names(rowDataSymbol) <- rowDataMaps[,1]
  rowDataSymbol <- rowDataSymbol[
    gsub("\\..*", "", rownames(bulkSummarizedExperiment))]
  length(rowDataSymbol)
  dim(bulkSummarizedExperiment)
  rowData(bulkSummarizedExperiment)$gene_symbol <- rowDataSymbol
  
  return(bulkSummarizedExperiment)
}


getExperimentData <- function(){
  # getExperimentData
  #
  # Get list of experiment data from paths.
  
  # get Y
  yPath <- "./data/monaco_et_al_2019/geo/GSE107011_Processed_data_TPM.txt"
  yList <- processY(yPath)
  
  # get Z
  zPath <- "./data/monaco_et_al_2019/manuscript/abisseq_rnaseq_cell_types_references_k17.csv"
  zref <- processZ(zPath)
  
  # get P
  pPath <- "./data/monaco_et_al_2019/manuscript/flow_cytometry_true_cell_type_proportions_s13.csv"
  ptrue <- processP(pPath, zref)
  
  returnList <- list(
    y.table = yList[["yTable"]],
    y.se = yList[["ySe"]],
    z = zref,
    p.true = ptrue
  )
  return(returnList)
}