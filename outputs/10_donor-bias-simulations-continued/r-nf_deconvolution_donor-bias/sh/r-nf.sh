#!/usr/bin/env sh

# Author: Sean Maden
#
# Run a nextflowr workflow.
#
#

# params
# directories
data_dir=data
rscript_dir=rscript
# workflow table path
wt_fname=workflow-table_intra.csv
wt_fpath=$data_dir/$wt_fname
# rscript paths
write_param_script=r-nf_write-params.R
gather_script=r-nf_gather-results.R

# update param.config
echo "updating params.config..."
Rscript ./$rscript_dir/$write_param_script -f $wt_fpath

# run main workflow
echo "running workflow..."
nextflow run main.nf

# gather results into table
echo "gathering results table..."
Rscript ./$rscript_dir/$gather_script