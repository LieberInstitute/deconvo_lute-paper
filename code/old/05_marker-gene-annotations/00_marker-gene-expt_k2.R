#!/usr/bin/env R
#
# Get marker genes for cell types using mean ratio expression (log counts). 
# Uses DeconvoBuddies::get_mean_ratio2() function.
#

# devtools::install_github("https://github.com/LieberInstitute/DeconvoBuddies")
libv <- c("DeconvoBuddies", "SingleCellExperiment", "SummarizedExperiment",
          "scater")
sapply(libv, library, character.only = T)

#-------
# params
#-------
celltype.varname <- "cellType_broad_hc"

#-----
# load
#-----
proj.dpath <- "deconvo_method-paper"

# path to full singlecellexperiment
sce.fname <- "sce_DLPFC.Rdata"
sce.fpath <- file.path("DLPFC_snRNAseq/processed-data/sce", sce.fname)
sce <- get(load(sce.fpath)) # get full singlecellexperiment

# save filepath
save.fnstem <- "k2"
sef.save.fname <- paste0("sef_mr-markers-",save.fnstem,
                         "-from-sce_dlpfc-ro1.rda")
scef.save.fname <- paste0("scef_mr-markers-",save.fnstem,
                          "-from-sce_dlpfc-ro1.rda")
markers.save.fname <- paste0("mr-markers-output_",save.fnstem,
                             "-ctbroadhc_from-sce_dlpfc-ro1.rda")
# make output filepaths
sef.fpath <- file.path(proj.dpath, "outputs/05_marker-gene-annotations",
                       sef.save.fname)
markers.fpath <- file.path(proj.dpath, "outputs/05_marker-gene-annotations",
                           markers.save.fname)
scef.fpath <- file.path(proj.dpath, "outputs/05_marker-gene-annotations",
                        scef.save.fname)

#---------------------
# summarize cell types
#---------------------
# variable name for cell types
table(sce[[celltype.varname]])
# Astro EndoMural     Micro     Oligo       OPC     Excit     Inhib
# 3979      2157      1601     32051      1940     24809     11067

# make celltypes variable
sce[["celltype"]] <- ifelse(grepl("Excit|Inhib", sce[[celltype.varname]]),
                            "Neuron", "Non-neuron")
table(sce[["celltype"]])
# Neuron Non-neuron
# 35876      41728

# summaries by slide
cd <- colData(sce)
cd$donor <- gsub("_.*", "", cd$Sample)
table(cd$donor, cd[,celltype.varname]) # full celltype variable
# cells by donor
df <- as.data.frame(table(cd$donor, cd[,"celltype"]))
df <- data.frame(donor = df[c(1:length(unique(cd$donor))),1],
                 neuron = df[df[,2]=="Neuron",3],
                 non_neuron = df[df[,2]=="Non-neuron",3])
df$sum <- df$neuron + df$glial
df$diff_g_minus_n <- df$glial - df$neuron
df
# slides by donor
cdf <- cd[!duplicated(cd$Sample),]
df <- as.data.frame(table(cdf$donor))
colnames(df) <- c("donor", "slides")
dff <- do.call(rbind, lapply(df[,1], function(di){
  di <- cd[cd$donor == di,]
  slidev <- unique(di$Sample)
  cellv <- unlist(lapply(slidev, function(si){
    nrow(di[di$Sample==si,])
  }))
  neuronv <- unlist(lapply(slidev, function(si){
    nrow(di[di$Sample==si & di$celltype == "Neuron",])
  }))
  glialv <- unlist(lapply(slidev, function(si){
    nrow(di[di$Sample==si & di$celltype == "Non-neuron",])
  }))
  num.slide <- length(slidev)
  return(c(sum(cellv)/num.slide, 
           sum(neuronv)/num.slide, 
           sum(glialv)/num.slide))
}))
colnames(dff) <- c("cells_per_slide", "neurons_per_slide", "non_neurons_per_slide")
df <- cbind(df, dff)


df$neurons_per_slide <- unlist(lapply(unique(cd$Sample), function(si){
  nrow(cdf[cdf$Sample==si,])
}))
df$glial_per_slide <- unlist(lapply(unique(cd$Sample), function(si){
  nrow(cdf[cdf$Sample==si,])
}))

#------------
# get markers
#------------
# get marker genes
markers <- get_mean_ratio2(sce, 
                           cellType_col = "celltype",
                           assay_name = "logcounts", 
                           add_symbol = TRUE)
# save results
save(markers, file = markers.fpath)

#------------------------------------------
# get sef -- subset top 100 markers by type
#------------------------------------------
# num. markers returned by cell type
table(markers$cellType.target)
# Neuron Non-neuron
# 4245        284

# get top 100 markers by gene
ctv <- unique(markers$cellType.target)
ma <- as.data.frame(markers, stringsAsFactors = F)
ma[,2] <- as.character(ma[,2])
markerv <- unique(unlist(lapply(ctv, function(ki){
  mi <- as.data.frame(ma[ma$cellType.target == ki,])
  mi <- mi[order(mi$rank_ratio),]; mi[seq(100),1]
})))

# make new se
scef <- sce[markerv,]
sef <- SummarizedExperiment(assays = list(counts = as.matrix(counts(scef)),
                                          logcounts = as.matrix(logcounts(scef))),
                            colData = colData(scef), rowData = rowData(scef))
# save sef
save(sef, file = sef.fpath)

# save new sec
scef <- scuttle::logNormCounts(scef)
save(scef, file = scef.fpath)

#------------------------------------------
# get sef -- subset top 20 markers by type
#------------------------------------------
num.markers <- 20
# get top 100 markers by gene
ctv <- unique(markers$cellType.target)
ma <- as.data.frame(markers, stringsAsFactors = F)
ma[,2] <- as.character(ma[,2])
markerv <- unique(unlist(lapply(ctv, function(ki){
  mi <- as.data.frame(ma[ma$cellType.target == ki,])
  mi <- mi[order(mi$rank_ratio),]; mi[seq(num.markers),1]
})))
# make new se
scef <- sce[markerv,]
sef <- SummarizedExperiment(assays = list(counts = as.matrix(counts(scef)),
                                          logcounts = as.matrix(logcounts(scef))),
                            colData = colData(scef), rowData = rowData(scef))
# manage new sef path
sef.save.fname <- paste0("sef_mr-markers_",save.fnstem,"_",
                         num.markers,"-per-k_dlpfc-ro1.rda")
sef.fpath <- file.path(proj.dpath, "outputs/05_marker-gene-annotations",
                       sef.save.fname)
# save sef
save(sef, file = sef.fpath)
