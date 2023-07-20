#---------------
# format results
#---------------
df.k2 <- rbind(df.s.k2.shared, df.s.k2.within)
# append abs.error
df.k2$abs.error.neuron <- abs(df.k2$neuron-df.k2$true.neuron)
df.k2$abs.error.glial <- abs(df.k2$glial-df.k2$true.glial)