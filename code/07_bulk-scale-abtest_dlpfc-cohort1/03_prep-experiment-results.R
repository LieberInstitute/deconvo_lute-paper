#---------------
# format results
#---------------
df.k2 <- rbind(df.s.k2.shared.counts, df.s.k2.within.counts,
               df.s.k2.shared.rpkm, df.s.k2.within.rpkm)
# append abs.error
df.k2$abs.error.neuron <- abs(df.k2$neuron-df.k2$true.neuron)
df.k2$abs.error.glial <- abs(df.k2$glial-df.k2$true.glial)