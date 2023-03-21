source("deconvo_method-paper/code/_halo-outputs-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
halo.output.table <- get(load(halo.output.path))
halo.output.table <- halo.output.table %>% as.data.frame()

# get transformed variables
# quantile transform marker counts by subject
levels.vector <- halo.output.table$Samples
halo.output.table[,marker.quantile.variable] <- levels.variable %>% 
  quantile_transform() %>% as.numeric()

# get inverse transformation
marker.vector <- halo.output.table[,gene.marker.label]
halo.output.table[,transformed.marker.variable] <- marker.vector %>%
  inverse_maximum_difference_transformation()

# resave outputs table
save(halo.output.table, file = output.updated.path)
