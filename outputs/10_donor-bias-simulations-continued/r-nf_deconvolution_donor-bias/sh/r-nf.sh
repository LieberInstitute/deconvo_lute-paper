#!/usr/bin/env sh

# Author: Sean Maden
#
# Run an r-nf workflow for deconvolution.
# 
# Flags:
# -w : Path to the workflow table.
# -p : Path to the parameter write .R script.
# -d : Path to the data directory.
# -r : Path to the rscripts directory.
# -n : Maximum number of runs per batch.
# 
# Returns:
# Runs the workflow, saving a new results table with unique timestamped filename on
# run success.
#

#------------
# parse flags
#------------
# initialize flag variables
workflow_table_filename=
params_write_rscript_filename=r-nf_write-params.R
results_gather_rscript_filename=r-nf_gather-results.R
data_directory=data
rscript_directory=rscript
runs_per_batch=200

#---------------------
# parse provided flags
#---------------------
while getopts wpgdrn: name; do
    case "$name" in
        'w')    workflow_table_filename=$OPTARG;;
        'p')    params_write_rscript_filename=$OPTARG;;
		'g')    results_gather_rscript_filename=$OPTARG;;
		'd')	data_directory=$OPTARG;;
		'r')	rscript_directory=$OPTARG;;
		'n')	runs_per_batch=$OPTARG;;
        \?)     echo "Invalid option provided: -$OPTARG" >&2;;
    esac
done

#-----------------
# manage arguments
#-----------------
workflow_table_path=$data_directory/$workflow_table_filename
write_param_script_path=$rscript_directory/$params_write_rscript_filename
results_gather_script_path=$rscript_directory/$results_gather_rscript_filename

#------------
# check paths
#------------
if [ ! -d "$workflow_table_path" ]; then
    echo 'ERROR: workflow table not found at provided path'
    exit 1
fi

if [ ! -d "$write_param_script_path" ]; then
    echo 'ERROR: write parameters .R script not found at provided path'
    exit 1
fi

if [ ! -d "$results_gather_script_path" ]; then
    echo 'ERROR: gather results .R script not found at provided path'
    exit 1
fi

#-----------------
# check runs count
#-----------------
run_count=$(cat $workflow_table_path | wc -l)
run_count=$(echo "$run_count-1" | bc)
batch_count=$(echo "$run_count/$num_batch" | bc)
echo 'Found '$run_count' runs.'
echo 'Parsing in '$batch_count' batches.'

#------------
# handle runs
#------------

#------------------
# aggregate results
#------------------

# update param.config
echo "updating params.config..."
Rscript ./$rscript_dir/$write_param_script -f $wt_fpath

# run main workflow
echo "running workflow..."
nextflow run main.nf

# gather results into table
echo "gathering results table..."
Rscript ./$rscript_dir/$gather_script