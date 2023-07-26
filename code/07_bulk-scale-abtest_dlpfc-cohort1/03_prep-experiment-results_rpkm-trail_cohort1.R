#---------------
# format results
#---------------
df.k2 <- rbind(df.s.k2.shared.counts.counts, df.s.k2.shared.rpkm.counts)
df.k2 <- rbind(df.k2, df.s.k2.shared.counts.logcounts)
df.k2 <- rbind(df.k2, df.s.k2.shared.rpkm.logcounts)
# append abs.error
df.k2$abs.error.neuron <- abs(df.k2$neuron-df.k2$true.neuron)
df.k2$abs.error.glial <- abs(df.k2$glial-df.k2$true.glial)