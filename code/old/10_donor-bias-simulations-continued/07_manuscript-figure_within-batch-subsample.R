source("deconvo_method-paper/code/10_donor-bias-simulations-continued/00_parameters.R")
sapply(libv, library, character.only = T)

# analyze results table
files.list <- list.files(data.dpath)
filter.files <- grepl(results.table.withingroup.filename.stem, files.list)
results.table.filename <- files.list[filter.files]
results.table.path <- file.path(save.dpath, rnf.dname, "data", results.table.withingroup.filename)
rt <- read.csv(rt.fpath)

method.varname <- "deconvolution_method"
# iterate on plot variables
varv <- c("prop.pred.type1", "prop.pred.type2", "bias.type1", 
          "bias.type2", "rmse.types")

# jitter plots
# filter duplicated data
rtf <- rt
rtf$label <- paste0(rtf$iterations_index, ";", rtf$deconvolution_method)
rtf <- rtf[!duplicated(rtf$label),]
methodv <- unique(rtf[,method.varname])
# get plots list
lgg <- lapply(varv, function(varname){
  dfp <- rtf; dfp$value <- dfp[,varname]
  ggplot(dfp, aes(x = deconvolution_method, y = value)) +
    geom_jitter(alpha = 0.5) + ggtitle(varname) +
    geom_boxplot(alpha = 0, color = "cyan", size = 1)
})
# bias plots
plot1 <- lgg[[3]] + theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  geom_hline(yintercept = 0, alpha = 1, color = "red")
plot2 <- lgg[[4]] + theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, alpha = 1, color = "red")
jpeg("ggjitterbox-comp-bias_intra-sample-bias.jpg", 
     width = 6, height = 5, units = "in", res = 400)
grid.arrange(plot1, plot2, nrow = 2,
             layout_matrix = matrix(c(1,1,1,2,2,2,2), ncol = 1),
             bottom = "Deconvolution method")
dev.off()

# plot pairs
# varname <- "prop.pred.type1"
# filter duplicate entries
lgg <- lapply(varv, function(varname){
  dfp <- do.call(cbind, lapply(methodv, function(methodi){
    filt <- rtf[,method.varname]==methodi; rti <- rtf[filt,]
    rti <- rti[order(rti$iterations_index),]; rti[,varname]
  }))
  dfp <- as.data.frame(dfp); colnames(dfp) <- methodv
  ggpairs(dfp) + ggtitle(varname)
})
names(lgg) <- varv

lgg$prop.pred.type1
lgg$prop.pred.type2
lgg$bias.type1
lgg$bias.type2
lgg$rmse.types

# violin plots
#
lgg <- lapply(varv, function(varname){
  dfp <- rtf; dfp$value <- dfp[,varname]
  ggplot(dfp, aes(x = deconvolution_method, y = value)) +
    geom_violin(draw_quantiles = 0.5) + ggtitle(varname)
})
# bias plots
plot1 <- lgg[[3]] + theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  geom_hline(yintercept = 0, alpha = 0.2)
plot2 <- lgg[[4]] + theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, alpha = 0.2)
jpeg("ggvp-comp-bias_intra-sample-bias.jpg", 
     width = 6, height = 5, units = "in", res = 400)
grid.arrange(plot1, plot2, nrow = 2,
             layout_matrix = matrix(c(1,1,1,2,2,2,2), ncol = 1),
             bottom = "Deconvolution method")
dev.off()

# scatter plot bias
# no facet
ggpt <- ggplot(rtf, aes(x = bias.type1, y = bias.type2, 
                        group = deconvolution_method, 
                        color = deconvolution_method)) +
  geom_point(alpha = 0.2) + geom_abline(slope = 1, intercept = 0) +
  theme_bw() + 
  geom_hline(yintercept = 0, alpha = 0.5, color = "gray") +
  geom_vline(xintercept = 0, alpha = 0.5, color = "gray")
# facet
ggpt <- ggplot(rtf, aes(x = bias.type1, y = bias.type2)) +
  geom_point(alpha = 0.2) + geom_abline(slope = 1, intercept = 0) +
  theme_bw() + 
  geom_hline(yintercept = 0, alpha = 0.5, color = "gray") +
  geom_vline(xintercept = 0, alpha = 0.5, color = "gray")
ggpt <- ggpt + facet_wrap(~deconvolution_method)
jpeg('ggpt-facet-method_intra-sample-bias.jpg', 
     width = 6, height = 6, 
     units = "in", res = 400)
ggpt
dev.off()