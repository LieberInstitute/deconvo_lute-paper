#!/usr/bin/env R

#
# Defines the cell size scale factors. Uses the following sets:
#
# * manul : the manually-defined set, also used in pseudobulk experiments.
# 
# * null : the null set, where factors are equal across types (both as "1").
#

# get cell size factor series
list.s.pred <- list(s.set.manual = c("glial" = 3, "neuron" = 10),
                    s.set.null = c("glial" = 1, "neuron" = 1))