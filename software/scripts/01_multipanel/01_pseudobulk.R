#!/usr/bin/env R

# Author: Sean Maden
#
# Run this script from: ./deconvo_method-paper/
#
# Makes environment for multipanel plot of results from k2, k3, k4 experiments across sn cohorts.
#
#
#

#--------
# cohort1
#--------
# K2
# load
load("./cohort1/env/02_pseudobulk/01_k2.RData") # pseudobulk, cohort1
# assign variables
dfp.tall.c1.k2 <- dfp.tall
rm(dfp.tall)

# K3
# load
load("./cohort1/env/02_pseudobulk/02_k3.RData") # pseudobulk, cohort1
# assign variables
dfp.tall.c1.k3 <- dfp.tall
rm(dfp.tall)

# K4
# load
# load("./cohort1/env/02_pseudobulk/03_k4.RData") # pseudobulk, cohort1
# assign variables
# dfp.tall.c1.k4 <- dfp.tall
# rm(dfp.tall)

#--------
# cohort2
#--------
# K2
# load
load("./cohort2/env/01_pseudobulk/01_k2_mrb_script.RData") # pseudobulk, cohort2
# assign variables
dfp.tall.c2.k2 <- dfp.tall
rm(dfp.tall)

# K3
# load
load("./cohort2/env/01_pseudobulk/01_k3_mrb_script.RData") # pseudobulk, cohort2
# assign variables
dfp.tall.c2.k3 <- dfp.tall
rm(dfp.tall)

# K4
# load
# load("./cohort2/env/02_pseudobulk/01_k4_mrb_script.RData") # pseudobulk, cohort2
# assign variables
# dfp.tall.c2.k4 <- dfp.tall
# rm(dfp.tall)


#-----
# save
#-----
# save env
# purge env
envDiff <- setdiff(ls(), 
  c("dfp.tall.c1.k2", "dfp.tall.c1.k3", 
    "dfp.tall.c2.k2", "dfp.tall.c2.k3"))
rm(list=envDiff)
# save env
save.image(
  file = "software/env/01_multipanel/01_pseudobulk_script.RData")