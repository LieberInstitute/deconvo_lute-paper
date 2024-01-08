#!/usr/bin/env R

### Author: Sean Maden
###
### Get diagram of deconvolution classes from lute.
###
### 1. A. deconvolutionParam
### 2. B. referencebasedParam
### 3. C. independentbulkParam
### 4. 1. nnlsParam
### 5. 2. musicParam
### 6. 3. epicParam
### 7. 4. deconrnaseqParam
### 8. 5. music2Param
### 9. 6. music2Param
### 10. 7. bisqueParam
### 11. 8. scdcParam
###
###
###

libv <- c("DiagrammeR", "lute")
sapply(libv, library, character.only=T)

## load data table
nodeTableName <- "lute-deconvolution_transfer-learning-table.nodeTable"
tablePath <- 
  file.path("./data/lute-deconvolution_transfer-learning-table.nodeTable")
nodeTable <- read.table(tablePath)

###-----------
### parameters
###-----------
appendString <- "Param"
methodColumName <- "method_class"
parentClassColumnName <- "parent_classes"

###-----------------
### helper functions
###-----------------
addMethodEdges <- function(filterString, 
                           nodeTableName, 
                           edgesInput, 
                           edgesOutput,
                           filterType="in", 
                           endAppendName="Param",
                           inputNodeString="referencebasedParam",
                           filterColumn="parent_classes",
                           methodColumn="diagram"){
## adds edges for methods, where originating node is called nodeString
methodVector <- nodeTable[,methodColumn]
filter <- grepl(filterString, nodeTable[,filterColumn])
if(!filterType == "in"){filter <- !filter}
newEdgeOut <- unique(methodVector[filter])
## get formatted end node names/edges
## newEdgeOut <- tolower(unique(methodVector[filter]))
## newEdgeOut <- paste0(newEdgeOut, endAppendName) # append string
newEdgeIn <- rep(
  edgesinputNodeString, length(newEdgeOut))
return(
  list(
    edgesInput=c(edgesInput, newEdgeIn),
    edgesOutput=c(edgesOutput, newEdgeOut)))
}

###----------------
### get chart edges
###----------------
## parse nodeTable data
## initialize parent class nodes and edges

edgesInput <- c("deconParam", "deconParam", "referencebasedParam")
edgesOutput <- c(
  "referencebasedParam", "referencefreeParam", "independentbulkParam")

## parse final child class edges
## append independent bulk methods
chartsList <- addMethodEdges(filterString="independentbulkParam",
                           filterType="in", 
                           nodeTable=nodeTable,
                           edgesInput=edgesInput, 
                           edgesOutput=edgesOutput, 
                           edgesinputNodeString="independentbulkParam")
## append remaining methods
chartsList <- addMethodEdges(filterString="independentbulkParam", 
                           filterType="out", 
                           nodeTable=nodeTable,
                           chartsList[["edgesInput"]], 
                           chartsList[["edgesOutput"]], 
                           edgesinputNodeString="referencebasedParam")

## get edge strings
edgeString <- paste0(chartsList[["edgesInput"]], "->", 
                      chartsList[["edgesOutput"]], collapse=" ")

###----------------
### get chart nodes
###----------------
## get chart strings
nodeVector <- unique(
  c(chartsList[["edgesInput"]], chartsList[["edgesOutput"]]))
filter <- nodeVector %in% nodeTable[,methodColname]
## get node strings
## set nodes for parent classes
nodeVectorFilter <- nodeVector[!filter]
nodeStringParent <- paste0(
  "node [shape=box, color='black',","fontname=Courier]\n",
  paste0(nodeVectorFilter, collapse=";"),"\n")
## set nodes for methods
nodeVectorFilter <- nodeVector[filter]
nodeStringMethod <- paste0(
  "node [shape=oval, color='black', fillcolor='lightgray', ",
  "fontname=Courier, style=filled]\n",
  paste0(nodeVectorFilter, collapse=";"), "\n")
## final nodeString
nodeString <- paste0(
  "\n\n",nodeStringParent,"\n\n", nodeStringMethod)

###----------------------------------
### make chart from component strings
#----------------------------------
## get chart strings
chartString <- paste0("digraph new_flowchart {\n",
                       "graph [overlap=false, fonsize=10]",
                       nodeString, "\n", 
                       edgeString, "\n",
                       "}")

grViz(chartString)

message(chartString)


###-----------------------
### define chart functions
###-----------------------
addMethodEdges <- function(filterString, nodeTable, edgesInput, edgesOutput,
                             filterType="in", endAppendName="Param",
                             edgesinputNodeString="referencebasedParam",
                             filterColumn="parent_classes",
                             methodColumn="method_class"){
  ## adds edges for methods, where originating node is called nodeString
  methodVector <- nodeTable[,methodColumn]
  filter <- grepl(filterString, nodeTable[,filterColumn])
  if(!filterType == "in"){filter <- !filter}
  newEdgeOut <- unique(methodVector[filter])
  newEdgeIn <- rep(edgesinputNodeString, length(newEdgeOut))
  return(
    list(
      edgesInput=c(edgesInput, newEdgeIn),
      edgesOutput=c(edgesOutput, newEdgeOut)))
}

deconvolutionClassesChart <- function(
    nodeTableFilename=
      paste0("lute-deconvolution","_transfer-learning-table.nodeTable"),
    methodColname="method_class",
    parentClassColname="parent_classes",
    #edgesInput.node.start=c("deconvolutionParam", 
    #                     "deconvolutionParam", 
    #                     "referencebasedParam"),
    #edgesOutput.node.start=c("referencebasedParam", 
    #                      "referencefreeParam", 
    #                      "independentbulkParam"),
    edgesInputNodeStart=c("A", "A", "B"),
    edgesOutputNodeStart=c("B", "B", "C"),
    methodNodeFillColor="lightgray",
    methodNodeOutLineColor="black",
    methodNodeShape="oval",
    parentNodeFillColor="white",
    parentNodeOutColor="black",
    parentNodeShape="box",
    methodNodeFont="Courier",
    parentNodeFont="Courier"){
  require(DiagrammeR)
  ## load nodeTable
  nodeTableName <- 
    "lute-deconvolution_transfer-learning-table.nodeTable"
  path <- file.path(
    system.file("nodeTable", package="lute"), nodeTableName)
  nodeTable <- read.table(path)
  
  ## set chart edges
  edgesInput <- edgesInputNodeStart; edgesOutput <- edgesOutputNodeStart
  ## append independent bulk methods
  #chartsList <- addMethodEdges(filterString="independentbulkParam",
  #                           filterType="in", nodeTable=nodeTable,
  #                           edgesInput=edgesInput, edgesOutput=edgesOutput, 
  #                           edgesinputNodeString="independentbulkParam")
  chartsList <- addMethodEdges(filterString="C",
                                 filterType="in", 
                                 nodeTable=nodeTable,
                                 edgesInput=edgesInput, 
                                 edgesOutput=edgesOutput, 
                                 edgesinputNodeString="C")
  ## append remaining methods
  #chartsList <- addMethodEdges(filterString="independentbulkParam", 
  #                           filterType="out", nodeTable=nodeTable,
  #                           chartsList[["edgesInput"]], chartsList[["edgesOutput"]], 
  #                           edgesinputNodeString="referencebasedParam")
  chartsList <- addMethodEdges(filterString="C", 
                             filterType="out", 
                             nodeTable=nodeTable,
                             chartsList[["edgesInput"]], 
                             chartsList[["edgesOutput"]], 
                             edgesinputNodeString="B")
  ## set edge strings
  edgeString <- paste0(
    chartsList[["edgesInput"]], "->", chartsList[["edgesOutput"]], collapse=" ")
  
  ## set node labels and formats
  nodeVector <- unique(
    c(chartsList[["edgesInput"]], chartsList[["edgesOutput"]]))
  filter <- nodeVector %in% nodeTable[,methodColname]
  ## set nodes for parent classes
  nodeVectorFilter <- nodeVector[!filter]
  nodeStringParent <- paste0("node [shape=", parentNodeShape, ", ",
                               "color='", parentNodeOutlineColor, "',",
                               "fillcolor='", parentNodeFillcolor, "',",
                               "fontname=", parentNodeFont, ",",
                               "style=filled]\n",
                               paste0(nodeVectorFilter, collapse=";"), "\n")
  ## set nodes for methods
  nodeVectorFilter <- nodeVector[filter]
  nodeStringMethod <- paste0("node [shape=", methodNodeShape, ", ",
                               "color='", methodNodeOutlineColor,"', ",
                               "fillcolor='", methodNodeFillcolor,"', ",
                               "fontname=", methodNodeFont, ", ",
                               "style=filled]\n",
                               paste0(nodeVectorFilter, collapse=";"), "\n")
  ## get final nodeString
  nodeString <- paste0("\n\n",nodeStringParent,"\n\n", nodeStringMethod)
  
  ## make chart from component strings
  ## get chart strings
  chartString <- paste0("digraph new_flowchart {\n",
                         "graph [overlap=false, fonsize=10]",
                         nodeString, "\n", edgeString, "\n", "}")
  chart <- grViz(chartString)
  return(
    list(
      chart=chart, chartString=chartString))
}

## example
newChart <- deconvolutionClassesChart()
newChart$chart

##-----
## save
##-----
jpeg("deconvolutionParam_hierarchy_diagram.jpg",
     res = 400, units = "in", width = 4, height = 4)
newChart$chart
dev.off()
