
libv <- c("ggplot2", "ggforce", "gridExtra")
sapply(libv, library, character.only = TRUE)

set.seed(0) # seed for random colors

#----------
# load data
#----------
# ro1 dlpfc sce, markers
fname1 <- "list-scef_markers-k2-k3-k4_ro1-dlpfc.rda"
lscef1 <- get(load(fname1))

# mrb dlpfc sce, markers
fname2 <- "list-scef_markers-k2-k3-k4_mrb-dlpfc.rda"
lscef2 <- get(load(fname2))

#-----------------------------
# dispersion plots, counts_adj
#-----------------------------
for(proj in c("ro1", "mrb")){
  message(proj)
  if(proj == "ro1"){
    lscef <- lscef1
    dfp <- metadata(lscef[[1]])[["dispersion_fits"]][["dfp"]]
  } else{
    lscef <- lscef2
    dfp <- metadata(lscef[[1]])[["dispersion_fits"]][["dfp"]]
    dfp[dfp$celltype=="other",]$celltype <- "glial"
  }
  dfp <- dfp[dfp$assay == "counts_adj",]
  dfp$marker.type <- ifelse(grepl("bg", dfp$marker.type),"b.g.","marker")
  
  # boxplots
  # boxplots with zoom
  line.width <- 0.5
  fill.col <- "cyan"
  lgg <- lapply(unique(dfp$celltype), function(ci){
    dfpf <- dfp[dfp$celltype==ci,]
    plot <- ggplot(dfpf, aes(x = marker.type, y = disp)) + 
      geom_boxplot(fill = fill.col, lwd = line.width) +
      facet_zoom(ylim = c(0, 10)) + theme_bw() + ggtitle(ci) +
      theme(axis.title.y = element_blank(),
            axis.title.x = element_blank())
    if(!ci == "all"){
      plot <- plot + theme(axis.text.x = element_blank(), 
                           axis.title.x = element_blank())
    } else{
      plot <- plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
    plot
  })
  num <- 5
  lm <- matrix(c(rep(1, num), rep(2, num), rep(3, num+1)), ncol = 1)
  fname <- paste0("ggbox-dispersion-composite_k2_",proj,"-dlpfc.jpg")
  jpeg(file = fname, width = 3.5, height = 4, units = "in", res = 400)
  grid.arrange(lgg[[1]], lgg[[2]], lgg[[3]], ncol = 1, 
               layout_matrix = lm,
               bottom = "Gene set", left = "Estimated dispersion")
  dev.off()
}

#------------------------------
# jitter plots with summary bar
#------------------------------
pt.alpha <- 0.5
pt.size <- 2
lgg <- lapply(unique(dfp$celltype), function(ci){
  dfpf <- dfp[dfp$celltype==ci,]
  plot <- ggplot(dfpf, aes(x = marker.type, y = disp)) + 
    geom_jitter(alpha = pt.alpha, size = pt.size) +
    stat_summary()
    facet_zoom(ylim = c(0, 10)) + theme_bw() + ggtitle(ci) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank())
  if(!ci == "all"){
    plot <- plot + theme(axis.text.x = element_blank(), 
                         axis.title.x = element_blank())
  } else{
    plot <- plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  plot
})
num <- 5
lm <- matrix(c(rep(1, num), rep(2, num), rep(3, num+1)), ncol = 1)
fname <- "ggbox-dispersion-composite_k2_ro1-dlpfc.jpg"
jpeg(file = fname, width = 3.5, height = 4, units = "in", res = 400)
grid.arrange(lgg[[1]], lgg[[2]], lgg[[3]], ncol = 1, 
             layout_matrix = lm,
             bottom = "Gene set", left = "Estimated dispersion")
dev.off()

#--------------------
# mean-variance plots
#--------------------