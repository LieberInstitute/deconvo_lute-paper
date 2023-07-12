# dependencies
libv <- c("here", "lute", "dplyr", "ggplot2", "gridExtra", "SingleCellExperiment", 
          "SummarizedExperiment", "scran", "MultiAssayExperiment")
sapply(libv, library, character.only = TRUE)

# save path
prepped.data.path <- here("deconvo_method-paper", "outputs", "01_prepare-datasets")
save.path <- here("deconvo_method-paper", "outputs", "05_cell-count-versus-cell-size")

#-------------------
# dlpfc cohort1 data
#-------------------
# multi assay experiment path
mae.filename <- "mae_final.rda"
mae.path <- here("deconvo_method-paper", "outputs", "01_prepare-datasets", mae.filename)
# dlpfc markers path
sce.markers.list.path <- here("deconvo_method-paper", "outputs", "01_prepare-datasets", "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda")

#-------------------
# dlpfc cohort2 data
#-------------------

# helper functions
append_k_columns <- function(df.input, df.ct = NULL, celltype.variable = "cell_type"){
  df.ct <- data.frame(cell_type = c("Inhib", "Other", "Astro", "Endo", "Excit", "Oligo", "OPC", "Micro", "all"),
                      k2 = c("neuron", "NA", "glial", "glial", "neuron", "glial", "glial", "glial", "NA"),
                      k3 = c("Inhib", "NA", "glial", "glial", "Excit", "glial", "glial", "glial", "NA"),
                      k4 = c("Inhib", "NA", "non_oligo_glial", "non_oligo_glial", "Excit", "Oligo", 
                             "non_oligo_glial", "non_oligo_glial", "NA"))
  celltype.vector <- df.input[,celltype.variable]
  unique.cell.types <- unique(celltype.vector)
  df.ct.filter <- df.ct[df.ct$cell_type %in% unique.cell.types,]
  df.ct.new <- do.call(rbind, lapply(celltype.vector, function(cell.type.id){
    df.ct.filter[df.ct.filter[,"cell_type"]==cell.type.id,,drop = F]
  }))
  cbind(df.input, df.ct.new)
}

# cell type harmonization table
harmonize_celltype_tables_1to1 <- function(df.input1, df.input2, celltype.variable = "cell_type"){
  # harmonize when there is one unique celltype level per row in each input
  celltype.vector1 <- df.input1[,celltype.variable]
  unique.cell.types1 <- unique(celltype.vector1)
  #df2.cn.filter <- !colnames(df.input2)==celltype.variable
  #df2.cn.return <- colnames(df.input2)[df2.cn.filter]
  
  do.call(rbind, lapply(unique.cell.types1, function(cell.type.id){
    df1.filter <- df.input1[,celltype.variable]==cell.type.id
    df2.filter <- df.input2[,celltype.variable]==cell.type.id
    c(df.input1[df1.filter,], df.input2[df2.filter,])
  })) %>% as.data.frame()
}