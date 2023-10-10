
# get cell type proportions

se_cell_size <- function(se, celltype.variable, sample.id.variable, 
                         label = "", method = "cell_median"){
  library(dplyr)
  sample.id.vector <- unique(se[[sample.id.variable]])
  df.size <- do.call(rbind, lapply(sample.id.vector, function(sample.id){
    se.filter <- se[,se[[sample.id.variable]]==sample.id]
    se.counts <- assays(se.filter)[["counts"]]
    celltype.vector <- unique(se.filter[[celltype.variable]])
    df.iter <- do.call(cbind, lapply(celltype.vector, function(cell.type){
      median(colSums(se.counts[,se.filter[[celltype.variable]]==cell.type]))
    }))
    df.iter <- as.data.frame(df.iter)
    colnames(df.iter) <- celltype.vector
    df.iter
  })) %>% as.data.frame()
  df.size$sample.id <- rownames(df.size) <- sample.id.vector
  df.size$size.type.label <- label
  df.size$size.method <- method
  return(df.size)
}

df_cell_size <- function(dftype, celltype.variable, sample.id.variable, 
                         cellsize.variable = "cell_size", label = ""){
  library(dplyr)
  sample.id.vector <- unique(dftype[,sample.id.variable])
  dftype.size <- do.call(rbind, lapply(sample.id.vector, function(sample.id){
    dff <- dftype[dftype[,sample.id.variable]==sample.id,]
    dff.size <- as.numeric(dff[,cellsize.variable])
    dff.size.return <- as.data.frame(t(dff.size))
    colnames(dff.size.return) <- dff$cell_type
    dff.size.return
  })) %>% as.data.frame()
  dftype.size$sample.id <- rownames(dftype.size) <- sample.id.vector
  dftype.size$size.type.label <- label
  return(dftype.size)
}

list_dfp_wide_tall_size <- function(mae){
  # sn proportions
  sn <- mae[["snrnaseq.k2.all"]]
  sn.size <- se_cell_size(sn, label = "snrnaseq", "k2", "Sample")
  # rn proportions
  rn <- mae[["cell.sizes"]]
  rn <- rn %>% t() %>% as.data.frame()
  rn.size <- df_cell_size(rn[rn$k.label=="k2",], label = "rnascope", "cell_type", "sample_id")
  
  # bind tall plot data
  colnames.bind <- intersect(colnames(sn.size), colnames(rn.size))
  rn.size <- rn.size[,colnames.bind]
  sn.size <- sn.size[,colnames.bind]
  dfp.tall <- as.data.frame(rbind(sn.size, rn.size))
  # variables
  dfp.tall$br.region <- gsub(".*_", "", dfp.tall$sample.id)
  dfp.tall$subject.id <- gsub("_.*", "", dfp.tall$sample.id)
  
  # bind wide plot data
  colnames(sn.size) <- paste0("sn.", colnames(sn.size))
  colnames(rn.size) <- paste0("rn.", colnames(rn.size))
  intersect.sample.id <- intersect(rownames(sn.size), rownames(rn.size))
  sn.size <- sn.size[intersect.sample.id,]
  rn.size <- rn.size[intersect.sample.id,]
  dfp <- cbind(sn.size, rn.size)
  dfp$sample.id <- rownames(dfp)
  # variables
  dfp$br.region <- gsub(".*_", "", dfp$sample.id)
  dfp$subject.id <- gsub("_.*", "", dfp$sample.id)
  
  return(list(dfp.wide = dfp, dfp.tall = dfp.tall, 
              df.sn.size = sn.size, df.rn.size = rn.size))
}



