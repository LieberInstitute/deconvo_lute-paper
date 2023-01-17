#!/usr/bin/env R

# Author: Sean Maden
#
# Cell size bias adjustment experiments
#

libv <- "lute"
sapply(libv, library, character.only = TRUE)

#--------------------
# test using sd_denom
#--------------------
# define offset series
set.seed(0)
pdiffv <- seq(-2e-2, 2e-2, 2e-3)
ddiffv <- rnorm(n = 21, mean = 1e-2, sd = 1e-3)

lseries <- lapply(seq(pdiffv), function(ii){
  p1true <- 0.75 + pdiffv[ii]; p2true <- 1-p1true
  donor_marker_biasexpt(offsetv = offsetv[ii], 
                        P = c(p1true, p2true),
                        donor.adj.method = "sd_denom",
                        bounds_thresh = 50)
})

dfr <- do.call(rbind, lapply(lseries, function(ii){ii$dfres}))

# scatterplots of results
dfp <- dfr
dfp$celltype <- ifelse(dfp$type.index==1, "neuron", "non-neuron")
dfp$offset.type <- ifelse(dfp$offset < 10 , "low", "high")

ggpt <- ggplot(dfp, aes(x = prop.true, y = prop.pred, 
                color = celltype, shape = celltype)) + theme_bw() + 
  geom_point(alpha = 0.5) + 
  geom_abline(slope = 1, intercept = 0, color = "black") +
  theme()

ggpt + facet_wrap(~prop.type*offset.type)

#-----------------------------------------
# test using sctransform::correct_counts()
#-----------------------------------------
install.packages("sctransform")


