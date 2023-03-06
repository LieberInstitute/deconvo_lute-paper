libv <- c("lute", "SummarizedExperiment", "SingleCellExperiment", 
          "ggplot2", "gridExtra", "GGally")
sapply(libv, library, character.only = TRUE)

#--------------------
# load results tables
#--------------------
# get save path
code.dname <- "10_donor-bias-simulations-continued"
proj.dname <- "deconvo_method-paper"
save.dpath <- file.path(proj.dname, "outputs", code.dname)
# get base path
rnf.dname <- "r-nf_deconvolution_intra-sample-bias"
base.path <- file.path(save.dpath, rnf.dname, "data")
# read tables
stemv <- paste0(c("inter", "intra"), "-sample-bias")
base.stem <- "r-nf_deconvolution_"
dnv <- paste0(base.stem, stemv)
lrt <- lapply(dnv, function(dni){
  base.pathi <- file.path(save.dpath, dni, "data")
  fni <- list.files(base.pathi)
  fni <- fni[grepl("^results_table_.", fni)]
  path <- file.path(base.pathi, fni)
  read.csv(path)
})

#-----------
# make plots
#-----------
# get plot data
dfp1 <- lrt[[1]]
dfp1$experiment <- "between-batch"
dfp2 <- lrt[[2]]
dfp2$experiment <- "within-batch"
dfp3 <- rbind(dfp1, dfp2)
# make tall plot data
typev <- unlist(strsplit(dfp3$type_labels[1], ";"))
num.types <- seq(length(typev))
dfp <- do.call(rbind, lapply(num.types, function(ii){
  dfi <- data.frame(method = dfp3$deconvolution_method,
                    experiment = dfp3$experiment,
                    value = dfp3[,paste0("bias.type", ii)])
  dfi$type <- typev[ii]; dfi
}))

# make new plot objects
ggjitter <- ggplot(dfp, aes(x = method, y = value)) + theme_bw() +
  geom_hline(yintercept = 0, color = "red", size = 2) +
  geom_jitter(alpha = 0.5) + geom_boxplot(alpha = 0, color = "cyan", size = 1)
ggjitter + facet_wrap(~experiment+type)
