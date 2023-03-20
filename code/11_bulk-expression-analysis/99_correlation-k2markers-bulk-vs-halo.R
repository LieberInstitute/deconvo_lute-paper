# load data
source("deconvo_method-paper/code/11_bulk-expression-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
rse.markers <- get(load(rse.k2markers.filepath))
halo <- get(load(halo.outputs.path))


#-----------------------------------------------
# correlation of bulk and halo marker expression
#-----------------------------------------------
# prepare correlation matrix
dfcor <- data.frame(akt3_expr = dfp.med[dfp.med[,1]=="akt3_expr",4],
                    nuc_area = dfp.med[dfp.med[,1]=="nuc_area",4],
                    nuc_perim = dfp.med[dfp.med[,1]=="nuc_perim",4],
                    slide = dfp.med[dfp.med[,1]=="nuc_perim",3],
                    type = dfp.med[dfp.med[,1]=="nuc_perim",2])
# get plot object
lcor <- lapply(unique(dfcor$type), function(ti){
  mcor <- cor(dfcor[dfcor$type==ti,seq(3)], method = "spearman")
})

# make correlation heatmap
names(lcor) <- unique(dfcor$type)
lgg <- lapply(seq(3), function(ii){
  ggcorrplot(lcor[[ii]], type = "lower", lab = T, 
             title = names(lcor)[ii])
})