# gets vector of sample ids to filter based on rnascope confidence annotations
  
file.path <- "./data/09_quality/rnascope_quality_annotations.csv"

anno <- read.csv(file.path)

samples.vector.exclude <- c()

filter.star <- anno$Star=="Low"

samples.vector.exclude <- paste0(anno[filter.star,]$Sample, ";Star")

filter.circle <- anno$Circle=="Low"

samples.vector.exclude <- c(samples.vector.exclude, paste0(anno[filter.circle,]$Sample, ";Circle"))





