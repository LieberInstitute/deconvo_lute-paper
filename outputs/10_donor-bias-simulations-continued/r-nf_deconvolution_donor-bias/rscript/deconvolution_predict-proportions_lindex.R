#!/usr/bin/env R

# Author: Sean Maden
#
# Run deconvolution functions on an SCE object, from command line.
#

libv <- c("SingleCellExperiment", "SummarizedExperiment", "argparse")
sapply(libv, library, character.only = T)

#-----------------
# parse parameters
#-----------------
# rscript.dir <- "$launchDir/rscript"
# decon.utils.scriptfname = "deconvolution_utilities.R"
# decon.utils.rscriptpath <- file.path(rscript.dir, decon.utils.scriptfname)

#-----------------
# source utilities
#-----------------
# source(decon.utils.rscriptpath)

#------------------
# helper functions
#------------------
predict_proportions <- function(Z, Y, Ybias = NULL, S = NULL, 
                                strict.method = "nnls", method.args = "", 
                                verbose = FALSE, ...){
  # predict_proportions
  #
  # Z : signature matrix
  # Y : bulk matrix
  # S : cell scale factors
  # strict_method : deconvolution method
  # method.args : additional method arguments
  # verbose : whether to show verbose messages
  #
  if(verbose){message("Getting cell type proportion predictions...")}
  
  # format string arguments
  cond <- method.args %in% c("", "NA")|is.na(method.args)
  if(!cond){
    method.args <- paste0(",",method.args)
  } else{
    method.args <- ""
  }
  method <- tolower(strict.method)
  
  # parse methods
  if(method == "nnls"){
    require(nnls)
    if(verbose){message("Using method nnls...")}
    command.string <- paste0("nnls::nnls(Z, Y",method.args,")$x")
  
  } else if(method == "music"){
    require(MuSiC)
    if(is(S, "NULL")){stop("Error, need mexpr to run method: ", method)}
    if(verbose){message("Using method MuSiC...")}
    if(verbose){message("Setting variances by gene...")}
    Sigma <- matrix(0, ncol = 1, nrow = nrow(Z))
    method.args <- ",S = S, Sigma = Sigma, nu = 1e-10, iter.max = 100, eps = 0"
    command.string <- paste0("music.basic(Y = Y, X = Z",method.args,")$p.weight")
  
    } else if(method == "bisque"){
      require(BisqueRNA)
      subject.variable <- "Sample"
      run.mode <- "default"
      if(run.mode == "default"){
        # append data to meet minimum sample requirements
        # use overlaps
        # use original data
        cd <- colData(sce)
        cdf <- cdf2 <- as.data.frame(cd)
        # force order
        var.lvl <- unique(cdf[,celltype.variable])
        cdf[,celltype.variable] <- factor(cdf[,celltype.variable], 
                                          levels=var.lvl[order(var.lvl)])
        cdf2[,subject.variable] <- "rep"
        rownames(cdf2) <- paste0(rownames(cdf2), "_2")
        cdf <- rbind(cdf, cdf2)
        adf <- AnnotatedDataFrame(cdf)
        sc.mexpr <- sc.mexpr2 <- as.matrix(assays(sce)[[assay.name]])
        colnames(sc.mexpr2) <- paste0(colnames(sc.mexpr2), "_2")
        sc.mexpr <- cbind(sc.mexpr, sc.mexpr2)
        sc.eset <- ExpressionSet(assayData = sc.mexpr, phenoData = adf)
        bulk.mexpr <- cbind(ypb, ypb, ypb)
        colnames(bulk.mexpr) <- c(unique(cdf[,subject.variable]), "2")
        bulk.eset <- ExpressionSet(assayData = bulk.mexpr)
        command.str <- paste0("ReferenceBasedDecomposition(",
                              "bulk.eset, sc.eset,",
                              "cell.types = celltype.variable,",
                              "subject.names = subject.variable,",
                              "use.overlap = TRUE",
                              ")$bulk.props[,'2']")
      } else{
        stop("Didn't recognize run mode.")
      }
  } else{
    if(verbose){message("Returning unmodified point prediction outputs.")}
  }
  # get proportion predictions
  p <- eval(parse(text = command.string))
  if(sum(p) > 1){p <- p/sum(p)} # coerce point scale
  if(verbose){message("Completed proportion predictions.")}
  
  return(p)
}

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
  message("index matrix data successfully loaded.")
} else{
  stop("Error, index matrix data not found at ",mi.filepath)
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
t1 <- Sys.time()
if(deconvolution.method=='music'){

  # get cell size factors for music method
  S <- unlist(lapply(unique.types, function(typei){
    type.filter <- celltype.vector==typei
    mean(colSums(mexpr[,type.filter]))
  }))
  pred.proportions <- predict_proportions(Z = Z, Y = Y, S = S,
                                          strict.method = deconvolution.method,
                                          unique.celltypes = unique.types)
} else{
  pred.proportions <- predict_proportions(Z = Z, Y = Y, strict.method = deconvolution.method,)
}

time.run <- Sys.time() - t1
time.run <- as.numeric(time.run)
names(time.run) <- "time.run.sec"

#---------------
# return results
#---------------
results.vector <- c()

# append params
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

# append proportions
names(pred.proportions) <- paste0("prop.pred.type", 
                                  seq(length(pred.proportions)))
results.vector <- c(results.vector, pred.proportions)

# append time, timestamp
results.vector <- c(results.vector, time.run, "timestamp" = as.numeric(t1))

# coerce to matrix
results.table <- t(as.data.frame(results.vector, nrow = 1))

# save new results
ts <- as.character(as.numeric(t1)) # get timestamp
new.filename <- paste0("deconvolution-results_", as.numeric(t1), ".csv")
write.csv(results.table, file = new.filename, row.names = FALSE)
