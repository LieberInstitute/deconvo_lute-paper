

#
# revisits ancova for halo confounds, with the following revisions:
# * subsample cells such that residuals are minimal
# * take subsets of cell sizes by size ranges such that residuals are minimal
# * do pca to compare cells, samples
#
#

# dependencies
libv <- c("here", "lute", "dplyr", "ggplot2", "gridExtra", "SingleCellExperiment", "scuttle",
          "SummarizedExperiment", "scran", "DeconvoBuddies", "UpSetR", "DelayedArray")
sapply(libv, library, character.only = TRUE)

# save path
save.path <- here("deconvo_method-paper", "outputs", "20_snrnaseq-bulk-matched_training")

#--------------
# prepare halo data
#--------------
# image reference
halo.path <- "Human_DLPFC_Deconvolution/processed-data/03_HALO/halo_all.Rdata"
halo.all <- get(load(halo.path))

# filter other category
halo.all <- halo.all[!halo.all$cell_type == "Other",]

#------------------
# get ancova models
#------------------
model.data <- halo.all[,c("Nucleus_Area", "cell_type", "BrNum", "Position", "Round", "Combo")]
ancova1 <- anova(lm(Nucleus_Area ~. , data = model.data))

# subsample on area
summary(halo.all$Nucleus_Area)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 5.006  17.800  29.913  32.750  42.274 199.999

fract.res <- function(anova.result, res.index = 5){
  anova.result$`Sum Sq`[res.index]/sum(anova.result$`Sum Sq`)}

model.filter1 <- model.data[model.data$Nucleus_Area <= 30,]
ancova.filter1 <- anova(lm(Nucleus_Area ~. , data = model.filter1))
res.filter1 <- fract.res(ancova.filter1)

model.filter2 <- model.data[model.data$Nucleus_Area <= 50,]
ancova.filter2 <- anova(lm(Nucleus_Area ~. , data = model.filter1))
res.filter2 <- fract.res(ancova.filter2)

model.filter3 <- model.data[model.data$Nucleus_Area <= 50 & model.data$Nucleus_Area >= 30,]
ancova.filter3 <- anova(lm(Nucleus_Area ~. , data = model.filter3))
res.filter3 <- fract.res(ancova.filter3)

model.filter4 <- model.data[model.data$Nucleus_Area <= 100 & model.data$Nucleus_Area >= 50,]
ancova.filter4 <- anova(lm(Nucleus_Area ~. , data = model.filter4))
res.filter4 <- fract.res(ancova.filter4)

#-----------------------------
# sample ids with fewest cells
#-----------------------------
# fewest total cells

# fewest neurons

# fewest glial

#----------------------------------------------
# subsample cells of varying amounts, by sample
#----------------------------------------------

