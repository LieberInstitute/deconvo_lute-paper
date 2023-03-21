source("deconvo_method-paper/code/_halo-outputs-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
halo.outputs.table <- get(load(halo.output.path))
halo.outputs.table <- halo.outputs.table %>% as.data.frame()

# 
cell.area.vector <- halo.outputs.table[,cell.area.variable]
halo.outputs.table[,cell.area.log.variable] <- cell.area.variable %>% log10()

# quantile transform marker counts by subject
samples.vector <- halo.outputs.table$Sample
unique.samples <- samples.vector %>% unique()
transformed.counts <- unlist(lapply(unique.samples, function(si){
  message("working on sample: ", si)
  filter <- halo.outputs.table$Sample==si
  gene.marker.vector <- halo.outputs.table[filter, gene.marker.label]
  quantile.vector <-  gene.marker.vector %>% quantile(seq(0, 1, 1e-3))
  quantile.labels.vector <- as.numeric(gsub("%", "", names(quantile.vector)))
  quantile.labels.vector2 <- quantile.labels.vector[1:(length(quantile.labels.vector)-1)]
  new.vector <- rep("NA", length(gene.marker.vector))
  for(ii in seq(length(quantile.labels.vector2))){
    which.marker <- gene.marker.vector >= quantile.labels.vector[ii] & 
      gene.marker.vector < quantile.labels.vector[ii+1]
    new.vector[which.marker] <- quantile.labels.vector[ii+1]
  }
  return(new.vector)
}))
transformed.counts <- as.numeric(transformed.counts)
transformed.counts[is.na(transformed.counts)] <- median(transformed.counts, na.rm = T) # replace NAs
halo.outputs.table[,marker.copies.quantile.variable] <- transformed.counts %>% as.numeric()

# plot marker quantiles by subject
jpeg(halo.quantiles.jpg.file.name, width = 10, height = 5, units = "in", res = 400)
transformed.counts %>% density() %>% 
  plot(title = "Marker Quantiles Frequency", xlab = "Frequency")
dev.off()

# qnorm method 2
# library(preprocessCore)

dfh$akt3.qnorm2 <- as.numeric(normalize.quantiles(as.matrix(dfh$AKT3_Copies, ncol = 1)))

quantile_normalization <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){return(my_mean[my_index])}
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean = df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

#-----------------------------------------------------
# plot distributions, before and after transformations
#-----------------------------------------------------
dfp1 <- data.frame(untransformed = dfh$Nucleus_Area, transformed = dfh$nuc.area.log10,
                   variable = "Nucleus_Area")

dfp2 <- data.frame(untransformed = dfh$AKT3_Copies, transformed = dfh$akt3.copies.qscale,
                   variable = "AKT3_Copies")

dfp <- rbind(dfp1, dfp2)

for(i in seq(2)){dfp[,i] <- as.numeric(dfp[,i])}

p1 <- ggplot(dfp, aes(x = untransformed)) + geom_density() + theme_bw() + ggtitle("Untransformed") + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  facet_wrap(~variable)

p2 <- ggplot(dfp, aes(x = transformed)) + geom_density() + theme_bw() + ggtitle("Transformed") + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  facet_wrap(~variable)

fname <- "ggdens-composite_cellsize-halo_ro1-dlpfc.jpg"
jpeg(file.path(save.dpath, fname), width = 5, height = 3.5, units = "in", res = 400)
grid.arrange(p1, p2, nrow = 2, left = "Density", bottom = "Value")
dev.off()
