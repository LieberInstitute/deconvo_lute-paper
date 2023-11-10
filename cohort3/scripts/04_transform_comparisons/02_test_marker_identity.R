#!/usr/bin/env R

# Author: Sean Maden
#
# Use to test the marker label identity by inspecting duplicated legend color 
# labels for 17 cell types.
#

libv <- c("ComplexHeatmap")
sapply(libv, library, character.only = TRUE)
knitr::opts_chunk$set(echo = TRUE)
setwd("..")
setwd("..")
load("./env/04_transform/01_transform_script.RData")

source("./notebooks/04_transform/00_param.R")

heatmapTall <- scale(zTpm)
mapTable <- mapTableImmuneCells("abis17")
markerTable <- markerTableFromList(listMarkers)
listMarkersHarmonize <- markersHarmonize(heatmapTall, markerTable)
heatmapTall <- listMarkersHarmonize[["heatmapTall"]]
markerTable <- listMarkersHarmonize[["markerTable"]]

dim(heatmapTall)
head(heatmapTall)
head(markerTable)
intersect(unique(markerTable$cellType), colnames(heatmapTall))
length(intersect(unique(markerTable$cellType), colnames(heatmapTall)))
mapTable
head(rownames(heatmapTall))
head(markerTable$markerName)
length(intersect(rownames(heatmapTall), markerTable$markerName))
length(intersect(rownames(heatmapTall), markerTable$markerName))
identical(rownames(heatmapTall), markerTable$markerName)

annotationHeatmapList <- annotationHeatmapList(
  heatmapTall, mapTable, markerTable)

annotationHeatmapList$heatmap
