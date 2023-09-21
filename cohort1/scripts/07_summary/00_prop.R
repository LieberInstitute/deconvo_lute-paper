
# get cell type proportions

se_cell_prop <- function(se, celltype.variable, sample.id.variable, label = ""){
  library(dplyr)
  sample.id.vector <- unique(se[[sample.id.variable]])
  df.prop <- do.call(rbind, lapply(sample.id.vector, function(sample.id){
    se[,se[[sample.id.variable]]==sample.id][[celltype.variable]] %>%
      table() %>% prop.table() %>% as.matrix() %>% t() %>% as.data.frame()
  })) %>% as.data.frame()
  df.prop$sample.id <- rownames(df.prop) <- sample.id.vector
  df.prop$prop.type.label <- label
  return(df.prop)
}

df_cell_prop <- function(dftype, celltype.variable, sample.id.variable, 
                         cellcount.variable = "cell_count", label = ""){
  library(dplyr)
  sample.id.vector <- unique(dftype[,sample.id.variable])
  dftype.prop <- do.call(rbind, lapply(sample.id.vector, function(sample.id){
    dff <- dftype[dftype[,sample.id.variable]==sample.id,]
    dff.counts <- as.numeric(dff[,cellcount.variable])
    names(dff.counts) <- dff$cell_type
    dff.prop <- dff.counts/sum(dff.counts)
    as.data.frame(t(dff.prop))
  })) %>% as.data.frame()
  dftype.prop$sample.id <- rownames(dftype.prop) <- sample.id.vector
  dftype.prop$prop.type.label <- label
  return(dftype.prop)
}

list_dfp_wide_tall <- function(mae){
  # sn proportions
  sn <- mae[["snrnaseq.k2.all"]]
  sn.prop <- se_cell_prop(sn, label = "snrnaseq", "k2", "Sample")
  # rn proportions
  rn <- mae[["cell.sizes"]]
  rn <- rn %>% t() %>% as.data.frame()
  rn.prop <- df_cell_prop(rn[rn$k.label=="k2",], label = "rnascope", "cell_type", "sample_id")
  
  # bind tall plot data
  dfp.tall <- as.data.frame(rbind(sn.prop, rn.prop))
  # variables
  dfp.tall$br.region <- gsub(".*_", "", dfp.tall$sample.id)
  dfp.tall$subject.id <- gsub("_.*", "", dfp.tall$sample.id)
  
  # bind wide plot data
  colnames(sn.prop) <- paste0("sn.", colnames(sn.prop))
  colnames(rn.prop) <- paste0("rn.", colnames(rn.prop))
  intersect.sample.id <- intersect(rownames(sn.prop), rownames(rn.prop))
  sn.prop <- sn.prop[intersect.sample.id,]
  rn.prop <- rn.prop[intersect.sample.id,]
  dfp <- cbind(sn.prop, rn.prop)
  dfp$sample.id <- rownames(dfp)
  # variables
  dfp$br.region <- gsub(".*_", "", dfp$sample.id)
  dfp$subject.id <- gsub("_.*", "", dfp$sample.id)
  
  return(list(dfp.wide = dfp, dfp.tall = dfp.tall))
}



