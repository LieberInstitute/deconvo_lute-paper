#!/usr/bin/env R

# Author: Sean Maden
#

libv <- c("ggplot2")
sapply(libv, library, character.only = T)

#----------
# load data
#----------
# read output df objects
read.fname <- "df-merge-cellsize_halo-scecounts.rda"
read.dpath <- out.dpath <- file.path("deconvo_method-paper", "outputs", 
                                     "07_cell-size-estimates")
dfm <- get(load(file.path(read.dpath, read.fname)))

#------------
# format data
#------------
# convert to numeric
for(c in seq(36)){dfm[,c] <- as.numeric(dfm[,c])}

#----------------------
# analysis by cell type
#----------------------
ctv <- c("Excit", "Inhib", "Oligo")
lcor <- lapply(ctv, function(cti){
  dffm <- dfm[,grepl(cti, colnames(dfm))]
  for(c in seq(ncol(dffm))){dffm[,c] <- as.numeric(dffm[,c])}
  cor(dffm, use = "pairwise.complete.obs", method = "spearman")
})

#-----------------------------
# compare means across samples
#-----------------------------
ctv <- c("Excit", "Inhib", "Oligo")
meanv <- colMeans(dfm[,seq(36)], na.rm = T)
dfp <- do.call(rbind, lapply(ctv, function(cti){
  mfv <- meanv[grepl(cti, names(meanv))]
  mfv.halo <- mfv[grepl("^halo\\..*", names(mfv))]
  dfp <- as.data.frame(matrix(mfv.halo, ncol = 1))
  dfp$halo.var <- gsub(paste0("\\.", cti, "$"), "", names(mfv.halo))
  dfp$celltype <- cti
  dfp$sce.mean.across.samples <- mfv[grepl("^sce\\..*", names(mfv))]
  return(dfp)
}))
colnames(dfp)[1] <- "mean.across.samples"

# get ggplot
ggpt <- ggplot(dfp, aes(x = sce.mean.across.samples, 
                        y = mean.across.samples,
                        color = celltype)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, col = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggpt+facet_wrap(~halo.var)




