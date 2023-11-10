#!/usr/bin/env R

# Author: Sean Maden
#
# Make heatmaps with annotations on rows and columns for cell types.
#
#
#

annotationHeatmapList <- function(
    heatmapTall, mapTable, markerTable){
  
  # column labels from 1 lookup
  columnLabelVector <- sapply(colnames(heatmapTall), function(typeLabel){
    mapTable[mapTable[,1]==typeLabel,2]})
  # row labels from 2 lookups
  rowLabelVector <- sapply(rownames(heatmapTall), function(markerLabel){
    typeLabel <- markerTable[
      markerTable[,"markerName"]==markerLabel,"cellType"]
    filterMapTable <- mapTable[,"cellType"]==typeLabel
    mapTable[filterMapTable,"color"]
  })
  
  # make annotations
  topAnnotation<-HeatmapAnnotation(
    typeLabel=colnames(heatmapTall),
    col=list(typeLabel=columnLabelVector)
  )
  #leftAnnotation<-rowAnnotation(
  #  markerTypeLabel=rownames(heatmapTall),
  #  col=list(markerTypeLabel=rowLabelVector)
  #)
  
  # match heatmap and marker table
  markerTable <- markerTable[
    order(match(markerTable$marker, rownames(heatmapTall))),]
  identical(rownames(heatmapTall), markerTable$marker)
  
  leftAnnotation<-rowAnnotation(
    markerTypeLabel=markerTable[,"cellType"],
    col=list(markerTypeLabel=columnLabelVector)
  )
  
  return(
    list(
      heatmapTall=heatmapTall,
      mapTable=mapTable,
      markerTable=markerTable,
      topAnnotation=topAnnotation,
      leftAnnotation=leftAnnotation,
      heatmap=Heatmap(
        heatmapTall, 
        top_annotation=topAnnotation, 
        left_annotation=leftAnnotation
      )
    )
  )
}

mapTableImmuneCells <- function(whichTable="abis17"){
  if(whichTable=="abis17"){
    mapTable <- data.frame(
      cellType = colnames(heatmapTall),
      color=c(
        "forestgreen", # Monocytes.C
        "tan", # NK
        "dodgerblue", # T.CD8.Memory
        "gold4", # T.CD4.Naive
        "gold4", # T.CD8.Naive
        "brown", # B.Naive
        "yellow", # T.CD4.Memory
        "blue4", # MAIT
        "purple", # T.gd.Vd2
        "gold2", # Neutrophils.LD
        "purple", # T.gd.non.Vd2
        "darkorange1", # "Basophils.LD
        "forestgreen", # Monocytes.NC.I
        "darkgoldenrod4", # B.Memory
        "darkgreen", # mDCs
        "orange", # pDCs
        "gold" # Plasmablasts
      )
    )
  } else{
    stop("Error, table not found.")
  }
  return(mapTable)
}

