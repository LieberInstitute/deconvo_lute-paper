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
    markerTypeLabel=markerTable$type,
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
