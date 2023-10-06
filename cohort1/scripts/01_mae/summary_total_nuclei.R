unique.samples <- unique(sc$Sample)
mean(unlist(lapply(unique.samples, function(sample.id){
  ncol(sc[,sc$Sample==sample.id])
})))

median(unlist(lapply(unique.samples, function(sample.id){
  ncol(sc[,sc$Sample==sample.id])
})))

sd(unlist(lapply(unique.samples, function(sample.id){
  ncol(sc[,sc$Sample==sample.id])
})))