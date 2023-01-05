#!/usr/bin/env R

# Author: Sean Maden
#
#
#

#----------
# load data
#----------
# read output df objects
read.fnv <- c("dfcellsize_halo.rda", 
              "df-cellsize_donor-region_sce-logcounts.rda")
read.dpath <- out.dpath <- file.path("deconvo_method-paper", "outputs", 
                                     "07_cell-size-estimates")
ldf <- lapply(read.fnv, function(fni){
  get(load(file.path(read.dpath, fni)))})
names(ldf) <- read.fnv

# save filename stem
save.fnstem <- "halo-scelct"

#----------------
# merge halo, sce
#----------------
dfs <- ldf[grepl("sce", names(ldf))][[1]]
dfh <- ldf$dfcellsize_halo.rda
# harmonize region variable
dfh$region <- gsub("[0-9]", "", dfh$donor)
dfh$region <- ifelse(dfh$region=="A", "ant",
                     ifelse(dfh$region == "M", "mid",
                            ifelse(dfh$region == "P", "post", "NA")))
# harmonize id variable
dfh$donor <- paste0("Br", gsub("[a-zA-Z]", "", dfh$donor))
# make new sample var
dfh$sample <- paste0(dfh$donor, "_", dfh$region)
# get overlapping donors, regions
sampv <- unique(intersect(dfh$sample, dfs$sample))
# filter df
dffh <- dfh[dfh$sample %in% sampv,]
dffs <- dfs[dfs$sample %in% sampv,]

# get new merged df, dfm
ctv <- c("Excit", "Inhib", "Oligo", "Micro", "OPC", "Astro")
colname.stem <- c("halo.akt3", "halo.cellarea", 
                  "halo.cytoarea", "halo.nucarea", 
                  "halo.nucparam", "sce.meanlct")
dfm <- do.call(rbind, lapply(sampv, function(sampi){
  dih <- dffh[dffh$sample == sampi,]
  dis <- dffs[dffs$sample == sampi,]
  dfm <- do.call(cbind, lapply(ctv, function(cti){
    # assign var defaults
    halo.akt3 <- halo.cellarea <- halo.cytoarea <- halo.nucarea <- 
      halo.nucparam <- sce.meancount <- "NA"
    # parse halo
    if(cti %in% dih$type){
      dihf <- dih[dih$type==cti,]
      halo.akt3 <- dihf[1]; halo.cellarea <- dihf[2]
      halo.cytoarea <- dihf[3]; halo.nucarea <- dihf[4];
      halo.nucparam <- dihf[5]
    }
    # parse sce
    if(cti %in% dis$celltype){
      disf <- dis[dis$celltype==cti,]
      sce.meancount <- disf$`mean(total_count)`[1]
    }
    matrix(unlist(c(halo.akt3, halo.cellarea, halo.cytoarea, halo.nucarea, 
                          halo.nucparam, sce.meancount)), nrow = 1)
  }))
  dfm <- as.data.frame(dfm)
  colnames(dfm) <- c(paste0(colname.stem, ".", 
                            rep(ctv, each = length(colname.stem))))
  dfm$sample <- sampi
  dfm
}))

# save merged df
save.fname <- paste0("df-merge-cellsize_",save.fnstem,".rda")
save.fpath <- file.path(read.dpath, save.fname)
save(dfm, file = save.fpath)
