#!/usr/bin/env R

# Author: Sean Maden
#
#

# load
mrna.yield <- read.csv(
  "./data/02_abisseq/s13_sample_proportions_mrna_yield.csv")
fc.proportions <- read.csv(
  "./data/02_abisseq/s13_sample_proportions_flow_cytometry.csv")
load("./env/02_abisseq/01_abisseq_script.RData")

# format
for(c in seq(2,ncol(fc.proportions))){
  fc.proportions[,c] <- as.numeric(fc.proportions[,c])
}
for(c  in seq(5,6)){
  mrna.yield[,c] <- as.numeric(mrna.yield[,c])
}

# save
save.image(file = "./env/02_abisseq/02_proportions_s13_script.RData")