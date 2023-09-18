#---------------
# format results
#---------------
list.counts <- list(
  df.s.k2.shared.counts.counts[,c(1:6,8:11)],
  df.s.k2.shared.rpkm.counts[,c(1:6,8:11)],
  df.s.k2.within.counts.counts,
  df.s.k2.within.rpkm.counts
)

list.logcounts <- list(
  df.s.k2.shared.counts.logcounts[,c(1:6,8:11)],
  df.s.k2.shared.rpkm.logcounts[,c(1:6,8:11)],
  df.s.k2.within.counts.lognorm,
  df.s.k2.within.rpkm.lognorm
)

# append bulk conditions
cd <- colData(mae[["bulk.rnaseq"]])
list.counts <- lapply(list.counts, function(df){
  df$bulk.sample.id <- rownames(df)
  df$bulk.sample.condition <- cd[df$bulk.sample.id,]$expt_condition
  df$assay.name.lutearg <- "counts"
  df
})
list.logcounts <- lapply(list.logcounts, function(df){
  df$bulk.sample.id <- rownames(df)
  df$bulk.sample.condition <- cd[df$bulk.sample.id,]$expt_condition
  df$assay.name.lutearg <- "logcounts"
  df
})

# bind final df
df.k2 <- rbind(do.call(rbind, list.counts), do.call(rbind, list.logcounts))

# append abs.error
df.k2$abs.error.neuron <- abs(as.numeric(df.k2$neuron)-as.numeric(df.k2$true.neuron))
df.k2$abs.error.glial <- abs(as.numeric(df.k2$glial)-as.numeric(df.k2$true.glial))

# format 
df.k2$true.neuron <- as.numeric(df.k2$true.neuron)

# clear cache
rm(df.s.k2.shared.counts.counts)
rm(df.s.k2.shared.rpkm.counts)
rm(df.s.k2.within.counts.counts)
rm(df.s.k2.within.rpkm.counts)
rm(df.s.k2.shared.counts.logcounts)
rm(df.s.k2.shared.rpkm.logcounts)
rm(df.s.k2.within.counts.lognorm)
rm(df.s.k2.within.rpkm.lognorm)
rm(mae)
rm(mae.final)
rm(rse.counts)
rm(rse.rpkm)
rm(sce.iter)

# save
save(df.k2, file = "./outputs/05_bulk/k2_results.rda")