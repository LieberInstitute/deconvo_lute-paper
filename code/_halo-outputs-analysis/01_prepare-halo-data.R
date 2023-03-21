source("deconvo_method-paper/code/_halo-outputs-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
halo.output.table <- get(load(halo.output.path))
halo.output.table <- halo.output.table %>% as.data.frame()

# get transformed variables
# quantile transform marker counts by subject
halo.output.table[,marker.quantile.variable] <- halo.output.table$Samples %>% 
  quantile_transform() %>% as.numeric()

# get inverse transformation
halo.output.table[,transformed.marker.variable] <- halo.output.table[,gene.marker.label] %>%
  inverse_maximum_difference_transformation()

# resave outputs table
save(halo.output.table, file = output.updated.path)





