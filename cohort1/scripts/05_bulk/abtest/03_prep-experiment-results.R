#---------------
# format results
#---------------
df.k2.counts <- rbind(df.s.k2.shared.counts.counts[,c(1:6,8:11)],
                      df.s.k2.shared.rpkm.counts[,c(1:6,8:11)],
                      df.s.k2.within.counts.counts,
                      df.s.k2.within.rpkm.counts)
df.k2.counts$assay.name.lutearg <- "counts"

df.k2.logcounts <- rbind(df.s.k2.shared.counts.logcounts[,c(1:6,8:11)],
                      df.s.k2.shared.rpkm.logcounts[,c(1:6,8:11)],
                      df.s.k2.within.counts.lognorm,
                      df.s.k2.within.rpkm.lognorm)
df.k2.logcounts$assay.name.lutearg <- "logcounts"
df.k2 <- rbind(df.k2.counts, df.k2.logcounts)

# append abs.error
df.k2$abs.error.neuron <- abs(as.numeric(df.k2$neuron)-as.numeric(df.k2$true.neuron))
df.k2$abs.error.glial <- abs(as.numeric(df.k2$glial)-as.numeric(df.k2$true.glial))

# save
save(df.k2, file = "./outputs/05_bulk/k2_results.rda")