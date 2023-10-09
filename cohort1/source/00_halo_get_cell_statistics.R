
#--------------------------
# prepare cell counts table
#--------------------------
h <- halo.output.table
sample.id.variable <- "SAMPLE_ID"
sample.id.combo.variable <- "Sample"
cell.area.variable <- "Nucleus_Area"
batch.variable <- "Combo"
cell.scale.factor.function <- "median"

# set filters
star.filter <- h[,batch.variable]=="Star"
circle.filter <- h[,batch.variable]=="Circle"

# apply k labels
# k2
h$k2 <- h$k3 <- h$k4 <- "other"
h[h$cell_type %in% c("Excit", "Inhib"),]$k2 <- "neuron"
h[h$cell_type %in% c("Oligo", "Micro", "Endo"),]$k2 <- "glial"
# k3
h[h$cell_type %in% c("Excit"),]$k3 <- "Excit"
h[h$cell_type %in% c("Inhib"),]$k3 <- "Inhib"
h[h$cell_type %in% c("Oligo", "Micro", "Endo"),]$k3 <- "glial"
# k4
h[h$cell_type %in% c("Excit"),]$k4 <- "Excit"
h[h$cell_type %in% c("Inhib"),]$k4 <- "Inhib"
h[h$cell_type %in% c("Oligo"),]$k4 <- "Oligo"
h[h$cell_type %in% c("Micro", "Endo"),]$k4 <- "non_oligo_glial"

# get list of sample id vectors
list.sample.id.filter <- list(combined.filter = seq(nrow(h)),
                              star.filter = which(h[,batch.variable]=="Star"),
                              circle.filter = which(h[,batch.variable]=="Circle"))

# get cell scale factors
list.cell.sizes.halo <- lapply(list.sample.id.filter, function(filter.iter){
  h.iter <- h[filter.iter,]
  # k2
  df.k2 <- aggregate(h.iter[,c(cell.area.variable, sample.id.combo.variable, "k2")], 
                     by = list(sample.id = h.iter[,sample.id.combo.variable], k2 = h.iter[,"k2"]), 
                     FUN = cell.scale.factor.function)[,c(1:3)]
  # k3
  df.k3 <- aggregate(h.iter[,c(cell.area.variable, sample.id.combo.variable, "k3")], 
                     by = list(sample.id = h.iter[,sample.id.combo.variable], k3 = h.iter[,"k3"]), 
                     FUN = cell.scale.factor.function)[,c(1:3)]
  # k4
  df.k4 <- aggregate(h.iter[,c(cell.area.variable, sample.id.combo.variable, "k4")], 
                     by = list(sample.id = h.iter[,sample.id.combo.variable], k4 = h.iter[,"k4"]), 
                     FUN = cell.scale.factor.function)[,c(1:3)]
  list(k2 = df.k2, k3 = df.k3, k4 = df.k4)
})
names(list.cell.sizes.halo) <- c("combo", "star", "circle")
save.path.cell.size.halo <- here("deconvo_method-paper", "outputs", "01_prepare-datasets", 
                                 "list_halo-cell-sizes_dlpfc-cohort1.rda")
save(list.cell.sizes.halo, file = save.path.cell.size.halo)

# get cell amounts
list.cell.amounts.halo <- lapply(list.sample.id.filter, function(filter.iter){
  h.iter <- h[filter.iter,]
  # k2
  df.k2 <- table(h.iter$k2, h.iter[,sample.id.combo.variable]) %>% as.data.frame()
  colnames(df.k2) <- c("k2", "sample.id", "cells")
  # k3
  df.k3 <- table(h.iter$k3, h.iter[,sample.id.combo.variable]) %>% as.data.frame()
  colnames(df.k3) <- c("k3", "sample.id", "cells")
  # k4
  df.k4 <- table(h.iter$k4, h.iter[,sample.id.combo.variable]) %>% as.data.frame()
  colnames(df.k4) <- c("k4", "sample.id", "cells")
  list(k2 = df.k2, k3 = df.k3, k4 = df.k4)
})
names(list.cell.sizes.halo) <- c("combo", "star", "circle")
save.path.cell.size.halo <- here("deconvo_method-paper", "outputs", "01_prepare-datasets", 
                                 "list_halo-cell-amounts_k2-k3-k4_dlpfc-cohort1.rda")
save(list.cell.sizes.halo, file = save.path.cell.size.halo)

# get cell proportions
list.cell.amounts.halo <- lapply(list.sample.id.filter, function(filter.iter){
  h.iter <- h[filter.iter,]
  
  # k2
  df.k2 <- do.call(rbind, lapply(unique(h.iter[,sample.id.combo.variable]), function(sample.id){
    dfi <- table(h.iter[h.iter[,sample.id.combo.variable]==sample.id,]$k2) %>% prop.table() %>% as.data.frame()
    colnames(dfi) <- c("k2", "prop")
    dfi$sample.id <- sample.id
    dfi
  }))
  
  # k3
  df.k3 <- do.call(rbind, lapply(unique(h.iter[,sample.id.combo.variable]), function(sample.id){
    dfi <- table(h.iter[h.iter[,sample.id.combo.variable]==sample.id,]$k3) %>% prop.table() %>% as.data.frame()
    colnames(dfi) <- c("k3", "prop")
    dfi$sample.id <- sample.id
    dfi
  }))
  
  # k4
  df.k4 <- do.call(rbind, lapply(unique(h.iter[,sample.id.combo.variable]), function(sample.id){
    dfi <- table(h.iter[h.iter[,sample.id.combo.variable]==sample.id,]$k4) %>% prop.table() %>% as.data.frame()
    colnames(dfi) <- c("k4", "prop")
    dfi$sample.id <- sample.id
    dfi
  }))
  
  list(k2 = df.k2, k3 = df.k3, k4 = df.k4)
})
names(list.cell.sizes.halo) <- c("combo", "star", "circle")
save.path.cell.size.halo <- here("deconvo_method-paper", "outputs", "01_prepare-datasets", 
                                 "list_halo-cell-proportions_k2-k3-k4_dlpfc-cohort1.rda")
save(list.cell.sizes.halo, file = save.path.cell.size.halo)
