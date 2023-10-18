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

# get id mappings
map.vector1 <- gsub("\\.", " ", colnames(prop.scaled))
map.vector2 <- df.proportions$cell.type
common.id <- intersect(map.vector1, map.vector2)
df.map <- data.frame(p.true.id = common.id,
                     map.vector2 = common.id)

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
  df.tall$true.prop.sd.mrna.yield <- "NA"
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
}

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



#-----
# save
#-----
save.image("./env/02_abisseq/01_abisseq_script.RData")
