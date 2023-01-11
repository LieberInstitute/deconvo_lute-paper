#!/usr/bin/env R

# Author: Sean Maden
#
# Perform constrained simulations with lute using the following params:
# K = 2 types total
# S: s1 > s2 (with slight random variation)
# P: p1~0.75; p2~0.25 (with slight random variation)
#
#

libv <- c("lute", "SingleCellExperiment", "SummarizedExperiment", "ggplot2")
sapply(libv, library, character.only = T)

#----------
# load data 
#----------
save.dpath <- file.path("deconvo_method-paper", "outputs", 
                        "07_cell-size-estimates")

# se marker data, k2
sef.fname <- "sef_mr-markers_k2_20-per-k_dlpfc-ro1.rda"
sef.dpath <- file.path("deconvo_method-paper", "outputs", 
                       "05_marker-gene-annotations")
sef <- get(load(file.path(sef.dpath, sef.fname)))

#-------------------------
# simulations -- sce.csize
#-------------------------
set.seed(0)
# simulation params
perc.var <- 0.4 # perc sd for sizes across sims
num.sim <- 10 # total sims
# make ls, varying the sizes slightly
s1 <- 10 # size of neurons
s2 <- 2 # size of non-neurons
s1diff <- rnorm(n = 10, sd = perc.var*s1)
s2diff <- rnorm(n = 10, sd = perc.var*s2)
lsv <- lapply(seq(num.sim), function(ii){
  c(s1 + s1diff[ii], s2 + s2diff[ii])
})
# get signature matrix, z
sef[["donor"]] <- sef[["BrNum"]]
setf <- set_from_sce(sef, groupvar = "donor", 
                     method = "mean", assayname = "logcounts")
lct <- assays(setf)$logcounts
# make lgv
lgv <- lapply(seq(num.sim), function(ii){list(lct[,1], lct[,2])})
# make lpv
p1v <- seq(0.7, 0.8, 0.1/num.sim); p2v <- 1-p1v
lpv <- lapply(seq(num.sim), function(ii){c(p1v[ii], p2v[ii])})

# run simulations
lres <- decon_analysis(lgv = lgv, lpv = lpv, lsv = lsv)

#-----------------
# make scatterplot
#-----------------
# get plot data
dfr <- lres$dfres
dfr$prop_k1_pred <- dfr$bias1 + dfr$prop_k1
dfr$prop_k2_pred <- dfr$bias2 + dfr$prop_k2
dfp <- rbind(data.frame(prop_true = dfr$prop_k1,
                        prop_pred = dfr$prop_k1_pred,
                        expt_type = dfr$zs_transform,
                        celltype = rep("Neuron", nrow(dfr))),
             data.frame(prop_true = dfr$prop_k2,
                        prop_pred = dfr$prop_k2_pred,
                        expt_type = dfr$zs_transform,
                        celltype = rep("Non-neuron", nrow(dfr))))

# get new expt lvl labels
lvlstr.false <- "not scaling by cells f"
lvlstr.true <- "lute (with scaling by cells f)"

# get rmse to print
rmse.false <- sqrt(mean((dfp[dfp$expt_type==FALSE,]$prop_true-
                          dfp[dfp$expt_type==FALSE,]$prop_pred)^2))
rmse.true <- sqrt(mean((dfp[dfp$expt_type==TRUE,]$prop_true-
                          dfp[dfp$expt_type==TRUE,]$prop_pred)^2))
dfp$rmse <- ifelse(dfp$expt_type==TRUE, rmse.true, rmse.false)
df.rmse <- data.frame(expt_type = c(lvlstr.false, lvlstr.true),
                    rmse = c(format(rmse.false, digits = 2), 
                             format(rmse.true, digits = 2)))
df.rmse$xpos <- 0.25; df.rmse$ypos <- 0.95
df.rmse$hjustpos <- df.rmse$vjustpos <- 0
df.rmse$rmse <- paste0("RMSE: ", df.rmse$rmse)

# format expt_type variable
dfp$expt_type <- ifelse(dfp$expt_type == "TRUE", lvlstr.true, lvlstr.false)
dfp$expt_type <- factor(dfp$expt_type, levels = c(lvlstr.false, lvlstr.true))

#---------------------
# make new plot object
#---------------------
# new plot object
ggpt <- ggplot() + theme_bw() +
  geom_text(data = df.rmse, alpha = 0.8,
            mapping = aes(x = xpos, y = ypos, 
                          hjust = hjustpos, vjust = vjustpos,
                          label = rmse)) +
  geom_point(dfp, mapping = aes(x = prop_true, y = prop_pred, 
                      shape = celltype, color = celltype),
             alpha = 0.5, size = 3) + 
  geom_abline(intercept = 0, slope = 1) +
  xlim(0.2, 0.8) + ylim(0, 1) +
  scale_color_manual(labels = c("Neuron", "Non-neuron"), 
                     values = c("blue", "red")) +
  xlab("True cell composition (cc)") +
  ylab("Estimated cc")
  
ggpt + facet_grid(cols=vars(expt_type))

#---------------
# save new plots
#---------------
plot.fnstem <- "k2-n10-markers-perk"
# make new pdf
plot.fname <- paste0("ggpt-facet-byexpttype_", plot.fnstem, ".pdf")
pdf(file.path(save.dpath, plot.fname), width = 5.8, height = 2.4)
ggpt + facet_wrap(~expt_type)
dev.off()
# make new jpeg
plot.fname <- paste0("ggpt-facet-byexpttype_", plot.fnstem, ".jpg")
jpeg(file.path(save.dpath, plot.fname), width = 5.8, height = 2.4, 
     units = "in", res = 400)
ggpt + facet_wrap(expt_type~.)
dev.off()
