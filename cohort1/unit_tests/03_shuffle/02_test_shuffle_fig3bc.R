#!/usr/bin/env R
# Author: Sean Maden

load("./env/03_shuffle/01_figbc_script.RData")

# expected values
expected.input.samples.count <- 10
expected.output.samples.count <- expected.input.samples.count

#------------------------------------------------------
# test 1, number of samples passed to algorithm wrapper
#------------------------------------------------------

length(unique(sce.k2$Sample))==expected.input.samples.count

#--------------------------------------------------------
# test 2, number of samples output from algorithm wrapper
#--------------------------------------------------------

length(unique(dfp.tall$sample.id))==expected.output.samples.count