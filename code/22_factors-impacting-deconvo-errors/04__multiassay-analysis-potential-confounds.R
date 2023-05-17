#!/usr/bin/env R

# Author: Sean Maden
#
# Joint analysis of halo, bulk, and sn data, and potential confounds
#
#

#---------------------------------
# prepare processed confounds data
#---------------------------------
# load confounds
halo.confounds <- get(load("list-halo-confounds_dlpfc-train.rda"))
bulk.confounds <- get(load("list-bulk-confounds_dlpfc-train.rda"))
sn.confounds <- get(load("df-confounds-sn_dlpfc-train.rda"))

# load marker expression
bulk.marker.expr <- get(load("list-donormarker-bulkexpr_dlpfc-train.rda"))
sn.marker.expr <- get(load("list-donormarker-snexpr_dlpfc-train.rda"))
