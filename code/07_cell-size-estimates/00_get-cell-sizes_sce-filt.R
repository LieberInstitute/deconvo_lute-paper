#!/usr/bin/env R

# Get the cell sizes from an SCE object.
#

libv <- c("SingleCellExperiment", "SummarizedExperiment")
sapply(libv, library, character.only = T)

#----------
# load data
#----------
# load sce
sce.fname <- "sef_mr-markers-k7-from-sce_dlpfc-ro1.rda"
sce.fpath <- file.path("deconvo_method-paper", "outputs", 
                       "05_marker-gene-annotations", sce.fname)
sce <- get(load(sce.fpath))

# out filename stem
out.dpath <- file.path()

#----------------------
# sizes by donor/region
#----------------------
dvarname <- "Sample"
dv <- unique(sce[[dvarname]])

# make cell types
ctvarname <- "cellType_broad_hc"
ctv <- sef[[ctvarname]]

# get total counts overall
agg.all <- data.frame(total_count = colSums(assays(sef)$counts),
                  celltype = ctv) %>% group_by(celltype) %>%
  summarise(mean(total_count), .groups = "drop") %>% as.data.frame()

# total counts by donor/region
agg.dr <- do.call(rbind, lapply(dv, function(di){
  seff <- sef[,sef[[dvarname]]==di] # filter
  dfi <- data.frame(total_count = colSums(assays(seff)$counts),
                        celltype = seff[[ctvarname]]) %>% 
    group_by(celltype) %>% 
    summarise(mean(total_count), .groups = "drop") %>% 
    as.data.frame()
  dfi$sample <- di
  dfi
}))
agg.dr$donor <- gsub("\\_.*", "", agg.dr$sample)
agg.dr$region <- gsub(".*\\_", "", agg.dr$sample)

#----------------
# sizes by region
#----------------
cd <- colData(sef)
cd$region <- gsub(".*\\_", "", cd$Sample)
rv <- unique(cd$region)

agg.rgn <- do.call(rbind, lapply(rv, function(ri){
  seff <- sef[,cd$region==ri] # filter
  dfi <- data.frame(total_count = colSums(assays(seff)$counts),
                    celltype = seff[[ctvarname]]) %>% 
    group_by(celltype) %>% 
    summarise(mean(total_count), .groups = "drop") %>% 
    as.data.frame()
  dfi$region <- ri
  dfi
}))

# correlations
var1v <- c("mid", "ant", "post")
var2v <- c("ant", "post", "mid")
mcor.rgn <- do.call(rbind, lapply(seq(var1v), function(ii){
  # get groups
  grp1 <- agg.rgn[agg.rgn$region==var1v[ii],]
  grp2 <- agg.rgn[agg.rgn$region==var2v[ii],]
  grp1 <- grp1[order(match(grp1$region, grp2$region)),]
  # get correlations
  cti <- cor.test(grp1[,2], grp2[,2])
  data.frame(donor = di,
             region1 = var1v[ii],
             region2 = var2v[ii],
             celltypes = paste0(unique(grp1$celltype), collapse = ";"),
             est = round(cti$estimate, 2),
             punadj = format(cti$p.value, scientific = T, digits = 3))
}))


#----------------------------------
# sizes by donor, binned replicates
#----------------------------------
# get donors with replicates
cd <- colData(sef)
colData(sef)$region <- gsub(".*\\_", "", cd$Sample)
colData(sef)$donor <- gsub("_.*", "", cd$Sample)
dv <- unique(cd$Sample)
ddf <- unique(gsub("_.*", "", dv[duplicated(gsub("_.*", "", dv))] ))

# list agg by donor, agg region for each
agg.drep <- do.call(rbind, lapply(ddf, function(di){
  sei <- sef[,sef$donor==di]
  rv <- unique(sei$region)
  dfi <- do.call(rbind, lapply(rv, function(ri){
    seff <- sei[,sei$region==ri] # filter
    dfi <- data.frame(total_count = colSums(assays(seff)$counts),
                      celltype = seff[[ctvarname]]) %>% 
      group_by(celltype) %>% 
      summarise(mean(total_count), .groups = "drop") %>% 
      as.data.frame()
    dfi$region <- ri
    dfi
  }))
  dfi$donor <- di
  dfi
}))

# correlations by donor
mcor <- do.call(rbind, lapply(ddf, function(di){
  message(di)
  ai <- agg.drep[agg.drep$donor==di,]
  riv <- unique(ai$region)
  # get overlapping cell types
  region.ol <- as.data.frame(table(ai$celltype))
  region.ol <- region.ol[region.ol[,2]==2,1]
  ai <- ai[ai$celltype %in% region.ol,]
  # get groups
  grp1 <- ai[ai$region==riv[1],]
  grp2 <- ai[ai$region==riv[2],]
  grp1 <- grp1[order(match(grp1$celltype, grp2$celltype)),]
  # get correlations
  cti <- cor.test(grp1[,2], grp2[,2])
  data.frame(donor = di,
             region1 = riv[1],
             region2 = riv[2],
             celltypes = paste0(unique(ai$celltype), collapse = ";"),
             est = round(cti$estimate, 2),
             punadj = format(cti$p.value, scientific = T, digits = 3))
}))

