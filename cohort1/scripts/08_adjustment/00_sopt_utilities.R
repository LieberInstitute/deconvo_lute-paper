

# helper functions
# get series of s cell size factors (THIS SCRIPT, AND A FEW OTHERS)
dfs.series <- function(s.glial.series = seq(1, 20, 1)){
  s.neuron.series <- rev(s.glial.series)
  dfs.series <- do.call(rbind, lapply(seq(length(s.glial.series)), function(index1){
    do.call(rbind, lapply(seq(length(s.neuron.series)), function(index2){
      c("glial" = s.glial.series[index1], "neuron" = s.neuron.series[index2])
    }))
  })) %>% as.data.frame()
  #plot(dfs.series$glial, dfs.series$neuron)
  return(dfs.series)
}

get_df_true_list <- function(sce, sample.id.variable = "Sample", 
                             celltype.variable = "k2"){
  sample.id.vector <- unique(sce[[sample.id.variable]])
  list.df.true <- lapply(sample.id.vector, function(sample.id){
    filter.sce <- sce[[sample.id.variable]]==sample.id
    prop.true <- prop.table(table(sce[,filter.sce][[celltype.variable]]))
    prop.true <- as.data.frame(t(as.matrix(prop.true)))
    rownames(prop.true) <- "true_proportion"
    prop.true
  })
  names(list.df.true) <- sample.id.vector
  return(list.df.true)
}