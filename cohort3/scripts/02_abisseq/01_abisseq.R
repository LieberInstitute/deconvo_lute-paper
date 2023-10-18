#!/usr/bin/env R

# Author: Sean Maden
#
# Performs deconvolution with the ABIS-seq signature matrix reference.
#
#
#

libv <- c("lute", "dplyr")
sapply(libv, library, character.only = TRUE)

#-----
# load
#-----
zref <- read.csv("./data/zref_abisseq.csv")
ptrue <- read.csv("./data/ptrue_abis.csv")
load("./env/01_tpm_summaries/01_read_script.RData")

#--------------
# format inputs
#--------------
# zref
zref <- zref[!zref[,1]=="12:00 AM",]
rownames(zref) <- zref[,1]
zref <- zref[,c(2:ncol(zref))]

# se
rownames(se) <- rowData(se)$gene_symbol
filter.rows <- !(duplicated(rownames(se))|is.na(rownames(se)))
se <- se[filter.rows,]

# ptrue
df.proportions <- ptrue
df.proportions <- df.proportions[c(3:nrow(df.proportions)),]
colnames(df.proportions) <- c("cell.type", 
                              "flow.cyto.mean", 
                              "flow.cyto.sd",
                              "mrna.yield.mean",
                              "mrna.yield.sd")
for(c in c(2:5)){df.proportions[,c] <- as.numeric(df.proportions[,c])}

# format cell type names: "T gd non Vd2"
df.proportions[df.proportions$cell.type=="T gd non-Vd2",]$cell.type <- "T gd non Vd2"

# get collapsed cell type ids
list.cell.type.names <- list("T CD8 Memory", "Monocytes NC I")
list.cell.types.vectors <- list(
  c("T CD8 CM", "T CD8 EM"), c("Monocytes NC", "Monocytes I")
)
for(index in seq(list.cell.type.names)){
  new.cell.type.name <- list.cell.type.names[[index]]
  cell.types.vector <- list.cell.types.vectors[[index]]
  new.row <- apply(df.proportions[df.proportions$cell.type %in% cell.types.vector,c(2:5)],2,mean)
  new.row <- c(new.cell.type.name, new.row)
  df.proportions <- as.data.frame(rbind(df.proportions, new.row))
  for(c in c(2:5)){df.proportions[,c] <- as.numeric(df.proportions[,c])}
  # drop old ids
  df.proportions <- df.proportions[!df.proportions$cell.type %in% cell.types.vector,]
}


#-----------------------
# get cell scale factors
#-----------------------
# marker library size
s.vector <- colSums(zref)

#------------
# run deconvo
#------------

result.unscaled <- lute(
  z = as.matrix(zref), 
  y = as.matrix(assays(se)[["tpm"]]), 
  assay.name = 'tpm',
  typemarker.algorithm = NULL
)

result.scaled <- lute(
  z = as.matrix(zref), 
  y = as.matrix(assays(se)[["tpm"]]),
  s = s.vector,
  assay.name = 'tpm',
  typemarker.algorithm = NULL
)

#------------------
# make df.plot.tall
#------------------
prop.unscaled <- result.unscaled[[1]]@predictions.table
prop.scaled <- result.scaled[[1]]@predictions.table

#---------------------------------
# get common cell type id mappings
#---------------------------------
map.vector1 <- gsub("\\.", " ", colnames(prop.scaled))
map.vector2 <- df.proportions$cell.type
common.id <- intersect(map.vector1, map.vector2)
df.map <- data.frame(p.true.id = common.id,
                     map.vector2 = common.id)

# filter common ids
filter.columns <- gsub("\\.", " ", colnames(prop.scaled)) %in% df.map[,1]
prop.scaled <- prop.scaled[,filter.columns]
prop.unscaled <- prop.unscaled[,filter.columns]

# bind results
prop.unscaled$sample.id <- rownames(prop.unscaled)
prop.unscaled$type <- "unscaled"
prop.scaled$sample.id <- rownames(prop.scaled)
prop.scaled$type <- "scaled"
df.plot.tall.s13 <- as.data.frame(rbind(prop.scaled, prop.unscaled))

#------------------------------
# aggregate results proportions
#------------------------------
# get means on overlapping ids
df.plot.tall.mean <- df.plot.tall.s13 %>% 
  group_by(type) %>% summarise(across(everything(), mean)) %>%
  t() %>% as.data.frame()
colnames(df.plot.tall.mean) <- df.plot.tall.mean[1,]
df.plot.tall.mean <- df.plot.tall.mean[seq(2,nrow(df.plot.tall.mean)),]
head(df.plot.tall.mean)

# get sds
df.plot.tall.sd <- df.plot.tall.s13 %>% 
  group_by(type) %>% summarise(across(everything(), sd)) %>%
  t() %>% as.data.frame()
colnames(df.plot.tall.sd) <- df.plot.tall.sd[1,]
df.plot.tall.sd <- df.plot.tall.sd[seq(2,nrow(df.plot.tall.sd)),]
head(df.plot.tall.sd)

# combine
colnames(df.plot.tall.sd) <- paste0("sd.", colnames(df.plot.tall.sd))
colnames(df.plot.tall.mean) <- paste0("mean.", colnames(df.plot.tall.mean))
df.tall <- as.data.frame(cbind(df.plot.tall.mean, df.plot.tall.sd))
# append true
df.tall$true.prop.mean.flow.cyto <- 
  df.tall$true.prop.sd.flow.cyto <- 
  df.tall$true.prop.mean.mrna.yield <- 
  df.tall$true.prop.sd.mrna.yield <-
  df.tall$s.cell.size <- "NA"
for(type in df.map$p.true.id){
  filter.tall <- rownames(df.tall) == type
  df.tall[filter.tall,"true.prop.mean.flow.cyto"] <- 
    df.proportions[df.proportions[,1]==type,"flow.cyto.mean"]
  df.tall[filter.tall,"true.prop.sd.flow.cyto"] <- 
    df.proportions[df.proportions[,1]==type,"flow.cyto.sd"]
  df.tall[filter.tall,"true.prop.mean.mrna.yield"] <- 
    df.proportions[df.proportions[,1]==type,"mrna.yield.mean"]
  df.tall[filter.tall,"true.prop.sd.mrna.yield"] <- 
    df.proportions[df.proportions[,1]==type,"mrna.yield.sd"]
  df.tall[filter.tall,"s.cell.size"] <- s.vector[type]
}

# format columns as numeric
for(c in seq(ncol(df.tall))){df.tall[,c] <- as.numeric(df.tall[,c])}
df.tall$cell.type <- rownames(df.tall)
# format cell size
df.tall$s.cell.size.log <- log(df.tall$s.cell.size)

# get df.wide
df.wide <- rbind(
  data.frame(mean = df.tall$mean.unscaled, sd = df.tall$sd.unscaled, 
             type = rep("unscaled", nrow(df.tall)), cell.type = rownames(df.tall),
             true.prop.mean.flow.cyto = df.tall$true.prop.mean.flow.cyto,
             true.prop.mean.mrna.yield = df.tall$true.prop.mean.mrna.yield,
             true.prop.sd.flow.cyto = df.tall$true.prop.sd.flow.cyto,
             true.prop.sd.mrna.yield = df.tall$true.prop.sd.mrna.yield),
  data.frame(mean = df.tall$mean.scaled, sd = df.tall$sd.scaled, 
             type = rep("scaled", nrow(df.tall)), cell.type = rownames(df.tall),
             true.prop.mean.flow.cyto = df.tall$true.prop.mean.flow.cyto,
             true.prop.mean.mrna.yield = df.tall$true.prop.mean.mrna.yield,
             true.prop.sd.flow.cyto = df.tall$true.prop.sd.flow.cyto,
             true.prop.sd.mrna.yield = df.tall$true.prop.sd.mrna.yield)
)
df.wide <- as.data.frame(df.wide)
df.wide <- df.wide[!df.wide$cell.type=="sample.id",]

# format columns as numeric
for(c in c(1,2,5,6,7,8)){df.wide[,c] <- as.numeric(df.wide[,c])}

# append bias and error
df.wide$true.prop.mean.flow.cyto.format <- 
  df.wide$true.prop.mean.flow.cyto/sum(na.omit(df.wide$true.prop.mean.flow.cyto))
df.wide$bias.flow.cyto <- df.wide$true.prop.mean.flow.cyto.format-
  df.wide$mean
df.wide$error.flow.cyto <- abs(df.wide$bias.flow.cyto)



#-----
# save
#-----
save.image("./env/02_abisseq/01_abisseq_script.RData")
