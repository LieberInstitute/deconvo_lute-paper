#!/usr/bin/env R

# Get the cell sizes from HALO outputs.
#
#

library(ggplot2)

#----------
# load data
#----------
# manage paths
base.path <- "."
save.dname <- "07_cell-size-estimates"
save.dpath <- file.path(base.path, "deconvo_method-paper", "outputs", save.dname)

# load rnascope data
lcsv.fname <- "lcsv_halo.rda"
lcsv.fpath <- file.path(base.path, "deconvo_method-paper", 
                        "data", "halo", lcsv.fname)
if(!file.exists(lcsv.fpath)){
  fpath <- file.path(base.path, "Human_DLPFC_Deconvolution", "raw-data", "HALO")
  fnv <- list.files(fpath, pattern = ".*.csv$", recursive=T) # get long paths
  fnv <- fnv[!grepl("prelim", fnv)]; dirv <- c("CIRCLE", "STAR")
  lcsv <- lapply(dirv, function(diri){
    fnvi <- fnv[grepl(paste0(".*",diri,".*"), fnv)]
    lcsvi <- lapply(fnvi, function(fni){read.csv(file.path(fpath, fni))})
    names(lcsvi) <- gsub(".*\\/|\\.csv$", "", fnvi)
    return(lcsvi)
  })
  names(lcsv) <- dirv
  # save binary
  save(lcsv, file = lcsv.fpath)
} else{
  lcsv <- get(load(lcsv.fpath))
}

#-----------
# set params
#-----------
labels = c("Endo" = "CLDN5", "Astro" = "GFAP",
           "Inhib" = "GAD1", "Excit" = "SLC17A7",
           "Micro" = "TMEM119", "Oligo" = "OLIG2")

cnv.size <- c("cell_area" = "Cell.Area..??m..", 
              "cytoplasm_area" = "Cytoplasm.Area..??m..", 
              "nucleus_area" = "Nucleus.Area..??m..",
              "nucleus_perimenter" = "Nucleus.Perimeter..??m.", 
              "akt3_expr" = "AKT3..Opal.570..Copies")

# get flattened, filtered list
lc1 <- lapply(lcsv[[1]], function(ii){ii})
lc2 <- lapply(lcsv[[2]], function(ii){ii})
lc <- c(lc1[grepl("Final", names(lc1))], lc2[grepl("Final", names(lc2))])

#-------------------
# get size variables
#-------------------
# make a table with one row per cell
# retain only cell type and size columns, and sample/donor info
dfs <- do.call(rbind, lapply(seq(length(lc)), function(ii){
  message(ii)
  csvi <- lc[[ii]]; csv.fname <- names(lc)[ii]
  cnv <- colnames(csvi)
  # get cell type label
  cnv.type <- colnames(csvi)[grepl(".*Positive$", colnames(csvi))]
  cnv.type <- cnv.type[grepl(paste0(labels, collapse = "|"), cnv.type)]
  do.call(rbind, lapply(cnv.type, function(typei){
    message(typei)
    # get filter terms
    which.type <- which(csvi[,typei]==1) # celltype filter
    csvif <- csvi[which.type,]
    # get size variable
    new.row <- rep("NA", 8)
    names(new.row) <- c(names(cnv.size), "num.cells", "type", "csv.fname")
    for(cni in cnv.size){
      # message(cni)
      cni.filt <- grepl(as.character(cni), colnames(csvif))
      if(length(which(cni.filt))==1){
        datv <- as.numeric(csvif[,cni.filt])
        new.row[names(cnv.size[cnv.size==cni])] <- round(mean(datv),3)
      }
    }
    new.row["num.cells"] <- length(which.type)
    new.row["type"] <- names(labels[labels==gsub("\\..*", "", typei)])
    new.row["csv.fname"] <- csv.fname; new.row
  }))
}))
dfs <- as.data.frame(dfs)
for(c in seq(5)){dfs[,c] <- as.numeric(dfs[,c])}
dfs$slide <- unlist(lapply(dfs$csv.fname, function(fni){
  unlist(strsplit(fni, "_"))[3]
}))
dfs$position <- gsub("[0-9]", "", dfs$slide)
dfs$donor <- gsub("[A-Z]", "", dfs$slide)

# save 
save.fname <- "dfcellsize-byslide_halo.rda"
save(dfs, file = file.path(save.dpath, save.fname))

#-------------
# violin plots
#-------------
dfp <- dfsize
colnames(dfp) <- c("akt3", "cell_area", "cyto_area", "nuc_area", 
                   "nuc_peri", "type", "fname", "donor")
dfp$type <- as.factor(dfp$type)

lgg <- lapply(colnames(dfp)[1:5], function(cni){
  dfpi <- dfp[,c(cni, "type")]; dfpi$y = dfpi[,cni]
  ggvp <- ggplot(dfpi, aes(x = type, y = y, color = type)) + 
    geom_violin(draw_quantiles = 0.5) +
    theme_bw() + ylab(cni)
  ggvp
})
names(lgg) <- colnames(dfp)[1:5]
