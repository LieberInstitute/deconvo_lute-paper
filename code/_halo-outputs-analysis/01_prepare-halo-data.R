source("deconvo_method-paper/code/_halo-outputs-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
halo.output.table <- get(load(halo.output.path))
halo.output.table <- halo.output.table %>% as.data.frame()

# normalize marker counts
marker.vector <- halo.output.table[,gene.marker.label]
levels.vector <- halo.output.table[,levels.variable]

# normalization1
halo.output.table[,normalization.variable1] <- levels.variable %>% normalization1() %>% as.numeric()

# normalization2
halo.output.table[,normalization.variable2] <- marker.vector %>% normalization2() %>% as.numeric()

# resave outputs table
save(halo.output.table, file = output.updated.path)
