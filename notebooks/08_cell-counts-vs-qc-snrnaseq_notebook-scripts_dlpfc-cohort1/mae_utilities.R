
#
# MAE.day -- standards for MultiAssayExperiments and similar object class types
#

assay_by_object_list <- function(mae, name.vector, assay.names.catch = c("counts", "Nucleus_Area")){
  #
  # get_assay_name_list
  # 
  assay.name.list <- lapply(name.vector, function(name.iter){
    assays.list.iter <- names(assays(mae[[name.iter]]))
    assays.list.iter[assays.list.iter %in% assay.names.catch]
  })
  names(assay.name.list) <- name.vector 
  return(assay.name.list)
}

get_mae_column_quantity_bygroup <- function(mae, group.variable.name, ...){
  quantity.matrix <- do.call(rbind, lapply(group.name.list, function(group.name){
    mae.iter <- mae.iter
    get_mae_column_quantity(mae.iter, ...)
  }))
  quantity.df <- as.data.frame(quantity.matrix)
  
}

map_column_quantity_options <- function(object, object.label, column.variable.name){
  se.cond <- is(object, "SingleCellExperiment")|is(object, "SummarizedExperiment")|is(object, "ExpressionSet")
  if(se.cond){ # classes containing coldata
    quantity.df <- as.data.frame(table(colData(object)[,column.variable.name]))
    
  } else{
    stop("unrecognized object class.")
  }
  quantity.df$mae.object.name <- object.label
  return(quantity.df)
  return(FALSE)
}

get_mae_column_quantity <- function(mae, object.name.vector = c("sn1.rnaseq", "rnascope.image"), 
                                    column.variable.levels = c("glial", "neuron"),
                                    quantity.type = "counts", column.variable.name = "k2",
                                    return.type = "data.frame"){
  #
  # example:
  # get_mae_column_quantity(mae, c("sn1.rnaseq", "rnascope.image"), c("glial", "neuron"), "counts", "k2")
  # get_mae_column_quantity(mae, c("sn1.rnaseq", "rnascope.image"), c("glial", "neuron"), "proportions", "k2")
  #
  #
  #
  
  # check for valid object names
  object.name.vector <- object.name.vector[object.name.vector %in% names(mae)]
  column.quantity.list <- lapply(object.name.vector, function(object.name){
    quantity.df <- map_column_quantity_options(mae[[object.name]], object.name, column.variable.name)
    quantity.df$quantity.type <- quantity.type
    quantity.df$cell.type.variable.name <- column.variable.name
    return(quantity.df)
  })
  # parse quantity.df options
  if(!is(column.variable.levels, "NULL")){
    column.quantity.list <- lapply(column.quantity.list, function(quantity.df){
      quantity.df[quantity.df[,1] %in% column.variable.levels,]
    })
  }
  if(quantity.type=="proportions"){
    column.quantity.list <- lapply(column.quantity.list, function(quantity.df){
      quantity.df[,2] <- quantity.df[,2]/sum(quantity.df[,2])
      quantity.df
    })
  }
  # parse return options
  if(!return.type=="list"){
    return.object <- as.data.frame(do.call(rbind, column.quantity.list))
    colnames(return.object) <- c("cell.type", "quantity", "mae.object.name", 
                                 "quantity.type", "cell.type.variable.name")
  } else{
    return.object <- column.quantity.list
  }
  return(return.object)
}

quantitydf_to_platformcelltype <- function(quantity.df){
  # gets wide table of platform-celltypes from quantity.df
  platforms <- unique(quantity.df[,3])
  unlist(lapply(platforms, function(platform.iter){
    qdf.iter <- quantity.df[quantity.df[,3]==platform.iter,]
    new.row.iter <- qdf.iter[, 2]
    names(new.row.iter) <- paste0(platform.iter, ";", qdf.iter[,1])
    new.row.iter
  }))
}

get_mae_quantity_bygroup <- function(mae, group.variable.name = "sample.id", 
                                     object.name.vector = c("sn1.rnaseq", "rnascope.image"), 
                                     column.variable.levels = c("glial", "neuron"),
                                     quantity.type = "counts", column.variable.name = "k2",
                                     return.type = "data.frame"){
  #
  # 
  # get_mae_quantity_bygroup(mae, "sample.id")
  #
  #
  
  
  unique.groups <- unique(colData(mae)[,group.variable.name])
  
  list.qdf <- lapply(unique.groups, function(group.iter){
    filter.mae <- colData(mae)[,group.variable.name]==group.iter
    mae.iter <- mae[,filter.mae,]
    qdf.iter <- get_mae_column_quantity(mae.iter, 
                                        object.name.vector = object.name.vector, 
                                        column.variable.levels = column.variable.levels,
                                        quantity.type = quantity.type, 
                                        column.variable.name = column.variable.name,
                                        return.type = return.type)
    qdf.iter$group.id <- group.iter
    return(qdf.iter)
  })
  
  # bind for tall qdf
  qdf.all <- do.call(rbind, lapply(list.qdf, function(item){item}))
  qdf.all <- as.data.frame(qdf.all)
  
  # get wide df
  qdf.wide <- do.call(rbind, lapply(list.qdf, function(qdf.iter){
    qdf.wide.iter <- quantitydf_to_platformcelltype(qdf.iter)
    c(qdf.wide.iter, unique(qdf.iter$group.id))
  }))
  qdf.wide <- as.data.frame(qdf.wide)
  colnames(qdf.wide)[ncol(qdf.wide)] <- "group.id"
  
  # return
  list.return.all.qdf <- list(qdf.tall = qdf.all, qdf.wide = qdf.wide)
  return(list.return.all.qdf)
}