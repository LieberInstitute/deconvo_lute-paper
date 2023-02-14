#!/usr/bin/env R

#
#
#
#

libv <- c("SingleCellExperiment", "limma", "SummarizedExperiment", "scuttle",
          "ggplot2")
sapply(libv, library, character.only = TRUE)

#---------------
# manage paths
#---------------
# get save dpath
code.dname <- "09_manuscript"
proj.dname <- "deconvo_method-paper"
save.dpath <- file.path(proj.dname, "outputs", code.dname)

#--------------
# manage params
#--------------
markeri <- "k2"
handle.str <- "ro1-dlpfc"

#----------
# load data
#----------
sce.fname <- paste0("sce_marker-adj-",markeri,"_",handle.str,".rda")
sce.fpath <- file.path(save.dpath, sce.fname)
scei <- get(load(sce.fpath))

#-----------------------------
# means and variances by donor
#-----------------------------
#groupvar <- paste0(scei[["k2"]], ";", scei[["Sample"]])
#table(groupvar)
#scei[["groupvar"]] <- groupvar
#sce.var <- aggregateAcrossCells(scei, ids = scei[["groupvar"]], statistics = "var")

assayname <- "counts_adj"
varname <- "groupvar"
scei[[varname]] <- paste0(scei[["k2"]], ";", scei[["Sample"]])
unique.groups <- unique(scei[[varname]])
sce.var <- do.call(cbind, lapply(unique.groups, function(groupi){
  message("working on group level ", groupi)
  filt <- scei[[varname]] == groupi
  rowVars(assays(scei[,filt])[[assayname]])
}))

sce.mean <- do.call(cbind, lapply(unique.groups, function(groupi){
  message("working on group level ", groupi)
  filt <- scei[[varname]] == groupi
  rowMeans(assays(scei[,filt])[[assayname]])
}))

colnames(sce.mean) <- colnames(sce.var) <- unique.groups

# save
# variances
fname <- paste0("matrix-donor-type_variances_",markeri, "_", handle.str,".rda")
save(sce.var, file = file.path(save.dpath, fname))
# means
fname <- paste0("matrix-donor-type_means_",markeri, "_", handle.str,".rda")
save(sce.mean, file = file.path(save.dpath, fname))

#-------------------
# plot mean-variance
#-------------------
genev <- rownames(sce.mean)
groupv <- colnames(sce.mean)
dfp <- do.call(rbind, lapply(seq(ncol(sce.mean)), function(ii){
  dfpi <- data.frame(gene = genev, mean = sce.mean[,ii], var = sce.var[,ii])
  dfpi$group <- groupv[ii]; return(dfpi)
}))
dfp$celltype <- gsub(";.*", "", dfp$group)
dfp$donor <- gsub(".*;", "", dfp$group)

ggsm <- ggplot(dfp, aes(x = mean, y = var, color = donor, group = donor)) +
  geom_smooth() + facet_wrap(~celltype) + theme_bw() +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  scale_x_log10() + scale_y_log10()

# save jpg
save.dpath <- "."
fname <- "ggsmooth-composite_k2-types-by-donor_ro1-dlpfc.jpg"
jpeg(file = file.path(save.dpath, fname), width = 8, height = 4, 
     units = "in", res = 400)
print(ggsm); dev.off()

#---------------------------------
# get random background dispersion
#---------------------------------


