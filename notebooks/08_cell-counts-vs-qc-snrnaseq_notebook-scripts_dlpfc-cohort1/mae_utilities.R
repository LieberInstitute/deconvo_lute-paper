
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

plot_bias <- function(df.lm.quantity){
  column.index.vector <- seq(ncol(df.lm.quantity))
  df.res <- do.call(cbind, lapply(column.index.vector, function(col.index){
    df.lm.iter <- data.frame(return.variable = df.lm.quantity[,col.index])
    col.seq.include <- column.index.vector[column.index.vector==!col.index]
    df.lm.iter <- cbind(df.lm.iter, df.lm.quantity[,col.seq.include])
    lm(return.variable ~ . , data = df.lm.iter)$residuals
  }))
  df.res <- as.data.frame(df.res)
  colnames(df.res) <- colnames(df.lm.quantity)
  # plot residuals
  ggpairs(df.res)
}

new_mae_with_filters <- function(mae, platform.name = "rnascope.image",
                                 sample.id.variable = "Sample",
                                 assay.name = "Nucleus_Area",
                                 value.filter.vector = c(100, 75, 50),
                                 old.platform.names.keep = c("sn1.rnaseq",
                                                             "rnascope.image")){
  #
  # example:
  # new_mae_with_filters(mae)
  #
  #
  
  #value.filter.vector <- c(100, 75, 50)
  #platform.name <- "rnascope.image"
  #assay.name <- "Nucleus_Area"
  #sample.id.variable <- "Sample"
  #old.platform.names.keep <- c("sn1.rnaseq", "rnascope.image")
  new.platform.names <- paste0(platform.name, "_filter-", 
                               tolower(assay.name),"_value-", value.filter.vector)
  
  # list of old platform objects
  old.platform.object.keep <- lapply(old.platform.names.keep, function(platform.name){
    mae[[platform.name]]
  })
  names(old.platform.object.keep) <- old.platform.names.keep
  # list of new platform objects
  mae.data <- mae[[platform.name]]
  new.platform.object.list <- lapply(value.filter.vector, function(filter.value){
    #filter.value <- cell.size.filters[1]
    new.data.name <- paste0(platform.name, "_filter", filter.value)
    filter.data.columns <- assays(mae.data)[[assay.name]] >= filter.value
    filter.data.columns <- as.vector(filter.data.columns[1,])
    mae.data.filtered <- mae.data[,filter.data.columns]
    # mae[[new.data.name]] <- mae.data.filtered
    mae.data.filtered
  })
  names(new.platform.object.list) <- new.platform.names
  
  # get output experiment list
  # bind lists
  output.platform.list <- append(old.platform.object.keep, new.platform.object.list)
  names(output.platform.list) <- c(names(old.platform.object.keep),
                                   names(new.platform.object.list))
  experiment.list.out <- ExperimentList(output.platform.list)
  
  # get sample map
  samplemap.old.keep <- sampleMap(mae)
  filter.samplemap.old <- samplemap.old.keep$assay %in% old.platform.names.keep
  samplemap.old.keep <- samplemap.old.keep[filter.samplemap.old,]
  samplemap.new <- do.call(rbind, lapply(seq(length(new.platform.object.list)),
                                         function(index){
                                           name <- names(new.platform.object.list)[index]
                                           item <- new.platform.object.list[[index]]
                                           data.frame(assay = rep(name, ncol(item)),
                                                      colname = colnames(item),
                                                      primary = item[[sample.id.variable]])
                                           
                                           
                                         }))
  samplemap.out <- rbind(samplemap.old.keep, samplemap.new)
  samplemap.out <- DataFrame(samplemap.out)
  
  #  get coldata
  coldata.out <- colData(mae)
  
  # get output mae
  mae.out <- prepMultiAssay(ExperimentList = experiment.list.out, 
                            sampleMap = samplemap.out, colData = coldata.out)
  mae.out.final <- MultiAssayExperiment(mae.out$experiments, 
                                        mae.out$colData, mae.out$sampleMap)
  
  return(mae.out.final)
}


#
#
#
#
pca_results_format <- function(df){
  #
  # example:
  # pca_results_format(df)
  #
  library(ggplot2)
  library(GGally)
  matrix.pca <- as.matrix(df)
  results.pca <- prcomp(x = matrix.pca)
  summary.pca <- summary(pca)
  # get data to plot
  plot.data.pca <- as.data.frame(results.pca$x)
  colnames(plot.data.pca) <- paste0(colnames(plot.data.pca), 
                                    " (", round(summary.pca$importance[2,], 2), "%)")
  # plot object
  ggpairs.plot <- ggpairs(plot.data.pca)
  # return
  return.list <- list(input.df = df,
                      results.pca = results.pca,
                       summary.pca = summary.pca,
                       plot.data.pca = plot.data.pca,
                       ggpairs.plot = ggpairs.plot)
  return(return.list)
}

#
#
#
pca_scatterplots_bygroups <- function(pca.results.object){
  #
  #
  #
  #
  
  # get plot data and format
  # get pca results to plot
  dfp <- pca.results.quantity$plot.data.pca
  xaxis.lab <- colnames(dfp)[1]
  yaxis.lab <- colnames(dfp)[2]
  colnames(dfp)[1:2] <- c("pc1", "pc2")
  # get grouping variables
  df.input <- pca.results.quantity$input.df
  dfp$group <- colnames(df.input)
  dfp$platform <- gsub(";.*", "", dfp$group)
  dfp$celltype <- gsub(".*;", "", dfp$group)
  
  # get plot objects
  # plot first 2 pcs
  plot.group1 <- ggplot(dfp, aes(x = pc1, y = pc2, color = group)) + geom_point()
  plot.group2 <- ggplot(dfp, aes(x = pc1, y = pc2, color = platform)) + geom_point()
  plot.group3 <- ggplot(dfp, aes(x = pc1, y = pc2, color = celltype)) + geom_point()
  plot.group1 <- plot.group1 + xlab(xaxis.lab) + ylab(yaxis.lab)
  plot.group2 <- plot.group2 + xlab(xaxis.lab) + ylab(yaxis.lab)
  plot.group3 <- plot.group3 + xlab(xaxis.lab) + ylab(yaxis.lab)
  
  # return
  return.list <- list(plot1 = plot.group1, plot2 = plot.group2, plot3 = plot.group3)
  return(return.list)
}