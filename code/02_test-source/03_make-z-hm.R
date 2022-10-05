#!/usr/bin/env R

#
# Make some heatmaps from z datasets
# 

# library(pheatmap)
library(ComplexHeatmap)
library(gridtext)
library(RColorBrewer)

#------
# paths
#------
proj.dpath <- "deconvo_method-paper"
save.dpath <- file.path(proj.dpath, "outputs/02_test-source")
lz.fname <- "lz_mr_dlpfc-ro1.rda"

# heatmaps to save
#
# filenames key:
# hmch : heatmap from ComplexHeatmap
# withscale/noscale : whether scale() function used before plotting
# coldefault/colphm : whether color palette is ComplexHeatmap default, or 
#   palette similar to pheatmap
# ...vs... : rows and columns
# 
# heatmap filenames
# donor;type table (e.g. z summary)
hm.png.fname1.1 <- "hmch-mr_withscale_coldefault_genemarkers-vs-donortype_dlpfc-ro1.png"
hm.png.fname1.2 <- "hmch-mr_noscale_coldefault_genemarkers-vs-donortype_dlpfc-ro1.png"
hm.png.fname1.3 <- "hmch-mr_withscale_colphm_genemarkers-vs-donortype_dlpfc-ro1.png"
hm.png.fname1.4 <- "hmch-mr_noscale_colphm_genemarkers-vs-donortype_dlpfc-ro1.png"
# cell type table (e.g. zexpt final)
hm.png.fname2.1 <- "hmch-mr_withscale_coldefault_genemarkers-vs-celltype_dlpfc-ro1.png"
hm.png.fname2.2 <- "hmch-mr_noscale_coldefault_genemarkers-vs-celltype_dlpfc-ro1.png"
hm.png.fname2.3 <- "hmch-mr_withscale_colphm_genemarkers-vs-celltype_dlpfc-ro1.png"
hm.png.fname2.4 <- "hmch-mr_noscale_colphm_genemarkers-vs-celltype_dlpfc-ro1.png"

#-----
# load
#-----
lz <- get(load(file.path(save.dpath, lz.fname)))

#-----------------
# helper functions
#-----------------
new_hm <- function(hm.fpath, hmi, rowlab, collab, legend.name, rescale = T,
                   show.col.names = F, use.colpal = c("ch","phm"), seed.num = 2,
                   png.width = 5, png.height = 4, png.res = 400, png.units = "in",
                   save.fig = T, append.rowanno = F, rowanno = NA,
                   append.colanno = F, colanno = NA){
  # saves a new heatmap to disk, then returns hm object
  #
  #
  set.seed(seed.num) # set seed for random color palette
  if(rescale){hmi <- base::scale(hmi)}
  if(use.colpal == "ch"){
    hm <- Heatmap(hmi, row_title = rowlab, 
                  column_title = collab, 
                  row_labels = "none", 
                  show_row_names = F, 
                  show_column_names = show.col.names, 
                  name = legend.name, 
                  show_row_dend = F, 
                  show_column_dend = F)
  } else{
    pal <- RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")
    col.pal <- colorRampPalette(rev(pal))(100)
    hm <- Heatmap(hmi, row_title = rowlab, 
                  column_title = collab, 
                  row_labels = "none", 
                  show_row_names = F,
                  show_column_names = show.col.names,
                  name = legend.name, 
                  show_row_dend = F, 
                  show_column_dend = F, 
                  col = col.pal)
  }
  if(append.rowanno){hm <- hm + rowanno}
  if(append.colanno){hm <- hm + colanno}
  if(save.fig){
    png(hm.fpath, width = png.width, height = png.height, res = png.res, 
        units = png.units); print(hm); dev.off()
  }
  return(hm)
}

get_lhm <- function(){
  #-------------------
  # z.summary heatmaps
  #-------------------
  hmi <- t(lz$z.summary.filt)
  collab <- paste0("Marker genes (N = ", ncol(hmi), ")")
  rowlab <- paste0("Donor;type (N = ",nrow(hmi),")")
  # get row anno
  annov <- rownames(hmi)
  dfa.row <- data.frame(donor = gsub(";.*", "", annov), 
                    cell_type = gsub(".*;", "", annov))
  rownames(dfa.row) <- rownames(hmi)
  # get col anno
  lz1 <- as.data.frame(lz[[1]])
  lz1 <- lz1[order(match(lz1[,1], colnames(hmi))),]
  cond <- identical(colnames(hmi), lz1[,1])
  dfa.col <- data.frame(cell_type = lz1[,2])
  # "hmch-mr_withscale_coldefault_genemarkers-vs-donortype_dlpfc-ro1.png"
  hm1 <- new_hm(hm.fpath = file.path(save.dpath, hm.png.fname1.1), hmi = hmi, 
                rowlab = rowlab, collab = collab, 
                legend.name = "Scaled\nMR_dk", 
                rescale = T, show.col.names = F, use.colpal = "ch",
                append.rowanno = T, rowanno = rowAnnotation(df = dfa.row),
                append.colanno = T, colanno = columnAnnotation(df = dfa.col))
  # "hmch-mr_noscale_coldefault_genemarkers-vs-donortype_dlpfc-ro1.png"
  hm2 <- new_hm(hm.fpath = file.path(save.dpath, hm.png.fname1.2), hmi = hmi, 
                rowlab = rowlab, collab = collab, legend.name = "MR_dk",
                rescale = F, show.col.names = F, use.colpal = "ch",
                append.rowanno = T, rowanno = rowAnnotation(df = dfa.row),
                append.colanno = T, colanno = columnAnnotation(df = dfa.col))
  # "hmch-mr_withscale_colphm_genemarkers-vs-donortype_dlpfc-ro1.png"
  hm3 <- new_hm(hm.fpath = file.path(save.dpath, hm.png.fname1.3), hmi = hmi,
                rowlab = rowlab, collab = collab, legend.name = "Scaled\nMR_dk",
                rescale = T, show.col.names = F, use.colpal = "phm",
                append.rowanno = T, rowanno = rowAnnotation(df = dfa.row),
                append.colanno = T, colanno = columnAnnotation(df = dfa.col))
  # "hmch-mr_noscale_colphm_genemarkers-vs-donortype_dlpfc-ro1.png"
  hm4 <- new_hm(hm.fpath = file.path(save.dpath, hm.png.fname1.4), hmi = hmi,
                rowlab = rowlab, collab = collab, legend.name = "MR_dk",
                rescale = F, show.col.names = F, use.colpal = "phm",
                append.rowanno = T, rowanno = rowAnnotation(df = dfa.row),
                append.colanno = T, colanno = columnAnnotation(df = dfa.col))
  #-----------
  # z heatmaps
  #-----------
  hmi <- lz$z.final
  rowlab <- paste0("Marker genes (N = ", nrow(hmi), ")")
  collab <- paste0("Cell types (N = ",ncol(hmi),")")
  # "hmch-mr_withscale_coldefault_genemarkers-vs-celltype_dlpfc-ro1.png"
  hm5 <- new_hm(hm.fpath = file.path(save.dpath, hm.png.fname2.1), hmi = hmi,
                rowlab = rowlab, collab = collab, legend.name = "Scaled\nMR_k",
                rescale = T, use.colpal = "ch", show.col.names = T,
                png.width = 4)
  # "hmch-mr_noscale_coldefault_genemarkers-vs-celltype_dlpfc-ro1.png"
  hm6 <- new_hm(hm.fpath = file.path(save.dpath, hm.png.fname2.2), hmi = hmi,
                rowlab = rowlab, collab = collab, legend.name = "MR_k",
                rescale = F, use.colpal = "ch", show.col.names = T,
                png.width = 4)
  # "hmch-mr_withscale_colphm_genemarkers-vs-celltype_dlpfc-ro1.png"
  hm7 <- new_hm(hm.fpath = file.path(save.dpath, hm.png.fname2.3), hmi = hmi,
                rowlab = rowlab, collab = collab, legend.name = "Scaled\nMR_k",
                rescale = T, use.colpal = "phm", show.col.names = T, 
                png.width = 4)
  # "hmch-mr_noscale_colphm_genemarkers-vs-celltype_dlpfc-ro1.png"
  hm8 <- new_hm(hm.fpath = file.path(save.dpath, hm.png.fname2.4), hmi = hmi,
                rowlab = rowlab, collab = collab, legend.name = "MR_k",
                rescale = F, use.colpal = "phm", show.col.names = T,
                png.width = 4)
  return(list(hm1, hm2, hm3, hm4, hm5, hm6, hm7, hm8))
}

#---------------
# table previews
#---------------
head(lz$top.marker.data)
#   A tibble: 6 × 8
#   gene    celltype.target mean.target cellType   mean ratio rank_ratio anno_ratio         
#   <chr>   <fct>                 <dbl> <fct>     <dbl> <dbl>      <int> <chr>              
# 1 MYO3B   Inhib                 0.268 OPC      0.0137 19.7           1 Inhib/OPC: 19.65   
# 2 SLC27A6 Inhib                 0.287 OPC      0.0192 15.0           2 Inhib/OPC: 14.973  
# 3 LYPD6B  Inhib                 0.301 Excit    0.0230 13.1           3 Inhib/Excit: 13.098
# 4 BTBD11  Inhib                 0.882 OPC      0.0796 11.1           4 Inhib/OPC: 11.084  
# 5 KIT     Inhib                 0.479 OPC      0.0464 10.3           5 Inhib/OPC: 10.328  
# 6 RELN    Inhib                 0.506 Oligo    0.0535  9.47          6 Inhib/Oligo: 9.465 

lz$z.summary.filt[c(1:5),c(1:6)]
#          Br2720;Inhib Br2720;Oligo  Br2720;OPC Br2720;Excit Br2720;Astro Br2720;EndoMural
#PRDM16      0.05363985  0.026162791 0.035545024   0.02063185  0.714285714       0.03883495
#LINC01141   0.04980843  0.008098007 0.004739336   0.02901354  0.008658009       0.00000000
#CAMK2N1     4.37876300  2.756021595 0.251184834   3.27208253  0.722943723       1.07766990
#SLC2A1      0.08100712  0.067898671 0.054502370   0.06060606  0.086580087       1.29126214
#RNF220      0.81171319  2.121885382 0.244075829   0.82849774  0.298701299       1.11650485

head(lz$z.final) 
#                Inhib      Oligo        OPC      Excit      Astro   EndoMural      Micro
# PRDM16    0.03021573 0.06817406 0.10516230 0.05010537 1.66901799 0.100346347 0.07531230
# LINC01141 0.03996508 0.01019263 0.02095346 0.12377250 0.00996760 0.008435346 0.93945279
# CAMK2N1   3.17718411 1.65782348 0.70413447 6.14272354 0.73512287 0.892507249 0.60191519
# SLC2A1    0.09119722 0.06702696 0.07649495 0.12652314 0.08209031 1.416279121 0.03084584
# RNF220    0.67761535 6.65279170 0.69598082 1.88259509 0.76918172 1.597458214 0.41346299
# MROH7     0.24592304 0.07453613 0.35925900 0.75625324 1.63512882 0.066032289 0.06529336

#------------------
# save new heatmaps
#------------------
get_lhm()
