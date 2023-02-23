library(ggplot2)
library(gridExtra)
library(ggpubr)

dft <- as.data.frame(table(sce[[group.variable]]))

unique.types <- unique(sce[[celltype.variable]])
dft <- do.call(rbind, lapply(unique.types, function(typei){
  dfi <- as.data.frame(table(sce[,sce[[celltype.variable]]==typei][[group.variable]]))
  dfi$type <- typei; dfi
}))

levelv <- unlist(lapply(unique(dft$Var1), function(ti){
  sum(dft[dft[,1]==ti,2])
}))
dft$Var1 <- factor(dft$Var1, levels = unique(dft$Var1)[rev(order(levelv))])

# get plot objects
ggbp1 <- ggplot(dft, aes(x = Var1, y = Freq, fill = type)) + 
  geom_bar(stat = "identity") + theme_bw() + ylab("Cell count") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())
ggbp.leg <- get_legend(ggbp1)
ggbp1 <- ggbp1 + theme(legend.position = "none")
ggbp2 <- ggplot(dft, aes(x = Var1, y = Freq, fill = type)) + 
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() + ylab("Fraction of cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(), 
        legend.position = "none")
  
lm <- matrix(c(rep(1, 5), rep(2, 5), 3), nrow = 1)

  grid.arrange(ggbp1, ggbp2, ggbp.leg, nrow = 1,
               layout_matrix = lm, bottom = "Sample ID")
  