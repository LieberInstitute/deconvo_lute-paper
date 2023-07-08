source("deconvo_method-paper/code/_halo-outputs-analysis/00_parameters.R")
sapply(libv, library, character.only = T)
halo.outputs.table <- get(load(halo.output.path))
halo.outputs.table <- halo.outputs.table %>% as.data.frame()

# 
cell.area.vector <- halo.outputs.table[,cell.area.variable]
halo.outputs.table[,cell.area.log.variable] <- cell.area.vector %>% log10()



# plot marker quantiles by subject
jpeg(halo.quantiles.jpg.file.name, width = 10, height = 5, units = "in", res = 400)
transformed.counts %>% density() %>% 
  plot(title = "Marker Quantiles Frequency", xlab = "Frequency")
dev.off()



#-----------------------------------------------------
# plot distributions, before and after transformations
#-----------------------------------------------------
dfp1 <- data.frame(untransformed = dfh$Nucleus_Area, transformed = dfh$nuc.area.log10, variable = "Nucleus_Area")

dfp2 <- data.frame(untransformed = dfh$AKT3_Copies, transformed = dfh$akt3.copies.qscale, variable = "AKT3_Copies")

dfp <- rbind(dfp1, dfp2)

for(i in seq(2)){dfp[,i] <- as.numeric(dfp[,i])}

p1 <- ggplot(dfp, aes(x = untransformed)) + geom_density() + theme_bw() + ggtitle("Untransformed") + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) + facet_wrap(~variable)

p2 <- ggplot(dfp, aes(x = transformed)) + geom_density() + theme_bw() + ggtitle("Transformed") + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) + facet_wrap(~variable)

fname <- "ggdens-composite_cellsize-halo_ro1-dlpfc.jpg"
jpeg(file.path(save.dpath, fname), width = 5, height = 3.5, units = "in", res = 400)
grid.arrange(p1, p2, nrow = 2, left = "Density", bottom = "Value")
dev.off()
