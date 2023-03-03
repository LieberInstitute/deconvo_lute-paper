#!/usr/bin/env R

# Author: Sean Maden
#
# Run deconvolution functions on an SCE object, from command line.
#

libv <- c("lute", "SingleCellExperiment", "SummarizedExperiment", "argparse")
sapply(libv, library, character.only = T)

#--------------
# manage parser
#--------------
parser <- ArgumentParser() # create parser object

# single-cell expression dataset
parser$add_argument("-r", "--sce_filepath", type="character", 
                    default="./data/sce-example.rda",
                    help = paste0("Expression reference dataset, ",
                                  "either an SCE or SE object."))

# bulk expression dataset
parser$add_argument("-b", "--bulk_filepath", type="character", default="NA",
                    help = "Bulk expression dataset. If NA, use means from sce")

# path to the iterations matrix
parser$add_argument("-l", "--list_index_filepath", type="character", default="NA",
                    help = "Filepath to list object containing iterations metadata.")

# iterations index
parser$add_argument("-i", "--iterations_index", type="character", default="1",
                    help = "Index of current run.")

# deconvolution parameters
parser$add_argument("-d", "--deconvolution_method", type="character", default="nnls",
                    help="Deconvolution method to use.")

# name of assay data in sce to use
parser$add_argument("-a", "--assay_name", type="character", default="counts",
                    help="Name of assay in reference expression data to use.")

# name of the celltype variable in sce coldata
parser$add_argument("-c", "--celltype_variable", type="character", default="celltype",
                    help="Reference ColData variable containing cell type labels.")

# get flags as variables
args <- parser$parse_args() # get parser object
# parse provided arguments
sce.filepath <- args$sce_filepath
bulk.filepath <- args$bulk_filepath
li.filepath <- args$list_index_filepath
index <- args$iterations_index
deconvolution.method <- args$deconvolution_method
assay.name <- args$assay_name
celltype.variable <- args$celltype_variable

#-----------------------
# load data, with checks
#-----------------------
# load single-cell expression data
if(file.exists(sce.filepath)){
  sce <- get(load(sce.filepath))
  message("Reference data successfully loaded.")
} else{
  stop("Error, reference data not found at ",sce.filepath)
}

# load index matrix
if(file.exists(li.filepath)){
  li <- get(load(li.filepath))
  message("index list data successfully loaded.")
} else{
  stop("Error, index list data not found at ",li.filepath)
}

#-----------------------
# parse iterations index
#-----------------------
# get celltype metadata
cd <- colData(sce)
if(celltype.variable %in% colnames(cd)){
  celltype.vector <- cd[,celltype.variable]
  unique.types <- unique(celltype.vector)
} else{
  stop("Error, didn't find ", celltype.variable," in reference colData.")
}

# subset types
index <- as.numeric(index)
lindex <- li[[index]] # get index
indexv <- lindex$vindex # sce filter index vector
sce <- sce[,indexv] # subset sce on indices
message("After index filter, retained ", ncol(sce), " cells.")

# overwrite new celltype data
cd <- colData(sce)
celltype.vector <- cd[,celltype.variable]
unique.types <- unique(celltype.vector)
unique.types <- unique.types[order(unique.types)]

#-------------------------------
# get deconvolution data objects
#-------------------------------
# get signature matrix
if(assay.name %in% names(assays(sce))){
  mexpr <- assays(sce)[[assay.name]]
  Z <- do.call(cbind, lapply(unique.types, function(ci){
    rowMeans(mexpr[,celltype.vector==ci])
  }))
  colnames(Z) <- unique.types
  S <- c("glial" = 3, "neuron" = 10)
  Z <- sweep(Z, 2, S, "*")
} else{
  stop("Error, assay not found in reference data.")
}

# load bulk expression data
# get bulk data
if(bulk.filepath == "NA"|is.na(bulk.filepath)|is(bulk.filepath, "NULL")){
  message("Making pseudobulk from reference data...")
  Y <- matrix(rowMeans(mexpr), ncol = 1)
} else if(file.exists(bulk.filepath)){
  Y <- get(load(bulk.filepath))
} else{
  stop("Error, bulk data not found at path: ", bulk.filepath)
}

#------------------
# run deconvolution
#------------------
arguments <- list()
if(deconvolution.method == "music"){
  arguments = list("S" = c("glial" = 3, "neuron" = 10))
}
deconvolution.results <- run_deconvolution(Z = Z, Y = Y, 
                                           method = deconvolution.method,
                                           arguments = arguments)

#---------------
# return results
#---------------
# begin new results table row
results.vector <- c()
results.vector["launch_dir"] <- getwd()
results.vector["sce_filepath"] <- basename(sce.filepath)
results.vector["bulk_filepath"] <- basename(bulk.filepath)
results.vector["li_filepath"] <- basename(li.filepath)
results.vector["iterations_index"] <- index
results.vector["num_cells"] <- ncol(sce)
results.vector["num_genes"] <- nrow(sce)
results.vector["deconvolution_method"] <- tolower(deconvolution.method)
results.vector["method_arguments"] <- "NA"
results.vector["assay_name"] <- assay.name
results.vector["type_labels"] <- paste0(unique.types, collapse = ";")
results.vector["celltype_variable"] <- celltype.variable

# get results values
# get timestamp
timestamp <- as.numeric(Sys.time())
# append proportions
propvalues <- deconvolution.results$predictions
names(propvalues) <- paste0("prop.pred.type",seq(length(propvalues)))
results.vector <- c(results.vector, propvalues)
# append time, timestamp
time.run <- as.numeric(deconvolution.results$benchmark[[1]][[1]])
time.run <- as.character(time.run)
results.vector <- c(results.vector, time.run, "timestamp" = timestamp)

# coerce to matrix
results.table <- t(as.data.frame(results.vector, nrow = 1))
colnames(results.table) <- names(results.vector)

# save new results
new.filename <- paste0("deconvolution-results_", timestamp, ".csv")
write.csv(results.table, file = new.filename, row.names = FALSE)
