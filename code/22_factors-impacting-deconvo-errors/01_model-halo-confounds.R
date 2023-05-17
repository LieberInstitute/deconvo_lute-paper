

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
total.cells <- as.data.frame(table(halo.all$SAMPLE_ID))
head(total.cells[order(total.cells[,2]),])
# Var1  Freq
# 7  Br3942M_CIRCLE  7585
# 2  Br2720P_CIRCLE  9303
# 35 Br8667M_CIRCLE 11969
# 18 Br6471A_CIRCLE 12612
# 1  Br2720M_CIRCLE 12688
# 14   Br6432A_STAR 13210

# fewest neurons
neuron.labels.vector <- c("Excit", "Inhib")
halo.filter <- halo.all$cell_type %in% neuron.labels.vector
total.cells <- as.data.frame(table(halo.all[halo.filter,]$SAMPLE_ID))
head(total.cells[order(total.cells[,2]),])
# Var1 Freq
# 22 Br6522M_CIRCLE  574
# 14   Br6432A_STAR 1624
# 4  Br2743A_CIRCLE 2993
# 33 Br8667A_CIRCLE 3181
# 2  Br2720P_CIRCLE 3495
# 24 Br6522P_CIRCLE 3570

# fewest glial
halo.filter <- !halo.all$cell_type %in% neuron.labels.vector
total.cells <- as.data.frame(table(halo.all[halo.filter,]$SAMPLE_ID))
head(total.cells[order(total.cells[,2]),])
# Var1 Freq
# 8    Br3942M_STAR 1416
# 36   Br8667M_STAR 1942
# 12   Br6423P_STAR 3010
# 7  Br3942M_CIRCLE 3449
# 23   Br6522M_STAR 5031
# 32   Br8492P_STAR 5764

#----------------------------------------------
# subsample cells of varying amounts, by sample
#----------------------------------------------

