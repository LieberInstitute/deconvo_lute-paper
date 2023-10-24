#!/usr/bin/env R

# Author: Sean Maden
#
# Append FC proportions and get RMSE by sample for S13 PBMC samples.
#

#-----
# load
#-----
load("./env/02_abisseq/01_abisseq_script.RData")
load("./env/02_abisseq/02_proportions_s13_script.Rdata")

df1 <- df.plot.tall.s13
df2 <- fc.proportions

#-------------------------------------------------------------
# get fc proportions as fractions of reference cells available
#-------------------------------------------------------------
# subset cells
vector.available.cells <- colnames(df1)[1:15]
df2 <- df2[,c(1,which(colnames(df2) %in% vector.available.cells))]
# convert to fractions
total.proportions <- rowSums(df2[,c(2:ncol(df2))])
df2[,c(2:ncol(df2))] <- t(
  apply(
    df2[,c(2:ncol(df2))], 1, function(ri){
      ri/sum(ri)}))
# check
rowSums(df2[,c(2:ncol(df2))])==1

#---------------------------------------
# append true proportions to predictions
#---------------------------------------
# differentiate cell column names
colnames(df1)[1:15] <- paste0(colnames(df1)[1:15],".nnls.pred")
colnames(df2)[2:ncol(df2)] <- paste0(colnames(df2)[2:ncol(df2)], ".fc.true")
# format sample ids
df1$sample.id.format <- gsub("_.*", "", df1$sample.id)
df2$sample.id.format <- gsub("_.*", "", df2$Sample.Name)
intersecting.sample.id <- intersect(df2$sample.id.format, df1$sample.id.format)
# merge
df3 <- merge(df1, df2, by = "sample.id.format")

# spot checks
cell.type.iter.df1 <- colnames(df1)[1]
cell.type.iter.df2 <- colnames(df2)[2]
sample.id.iter <- intersecting.sample.id[1]
df1[df1$sample.id.format==sample.id.iter,
    cell.type.iter.df1]
df2[df2$sample.id.format==sample.id.iter,
    cell.type.iter.df2]
df3[df3$sample.id.format==sample.id.iter,
    cell.type.iter.df1]
df3[df3$sample.id.format==sample.id.iter,
    cell.type.iter.df2]

#---------------
# calculate RMSE
#---------------
# get bias, error, rmse

df.rmse.scaled <- df.rmse.unscaled <- data.frame(cell.type = vector.available.cells)
df.rmse.scaled$rmse <- df.rmse.unscaled$rmse <- "NA"

for(cell.type in vector.available.cells){
  message(cell.type)
 vector.pred <- df3[,grepl(cell.type, colnames(df3)) & grepl("nnls.pred", colnames(df3))]
 vector.true <- df3[,grepl(cell.type, colnames(df3)) & grepl("fc.true", colnames(df3))]
 df3[,ncol(df3)+1] <- vector.pred-vector.true
 colnames(df3)[ncol(df3)] <- paste0(cell.type,".bias.pred.true")
 df3[,ncol(df3)+1] <- abs(vector.pred-vector.true)
 colnames(df3)[ncol(df3)] <- paste0(cell.type,".error.pred.true")
 # append rmse
 df.rmse.scaled[df.rmse.scaled==cell.type,]$rmse <- sqrt(
   mean(
     df3[
       df3$type=="scaled",paste0(cell.type,".error.pred.true")]^2))
 df.rmse.unscaled[df.rmse.unscaled==cell.type,]$rmse <- sqrt(
   mean(
     df3[
       df3$type=="unscaled", paste0(cell.type,".error.pred.true")]^2))
}

df.rmse.scaled$type <- "scaled"
df.rmse.unscaled$type <- "unscaled"

# df.rmse final tables
df.rmse.tall <- as.data.frame(rbind(df.rmse.scaled, df.rmse.unscaled))
df.rmse.wide <- data.frame(
  cell.type = df.rmse.scaled$cell.type,
  rmse.scaled = df.rmse.scaled$rmse,
  rmse.unscaled = df.rmse.unscaled$rmse
)

# append cell sizes
df.rmse.wide$s.cell.size <- df.rmse.tall$s.cell.size <- "NA"
for(cell.type in df.rmse.wide$cell.type){
  cell.size.iter <- df.tall[df.tall$cell.type==cell.type,]$s.cell.size[1]
  df.rmse.wide[df.rmse.wide$cell.type==cell.type,]$s.cell.size <- 
    df.rmse.tall[df.rmse.tall$cell.type==cell.type,]$s.cell.size <-
    cell.size.iter
}

#-----
# save
#-----
save.image("./env/02_abisseq/03_rmse_script.RData")
