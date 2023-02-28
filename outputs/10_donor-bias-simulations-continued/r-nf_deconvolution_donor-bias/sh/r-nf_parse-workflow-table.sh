#!/bin/bash

# Author: Sean Maden
#
# Run an r-nf workflow for deconvolution. Since large channel lists cannot presently be
# processed all at once, workflow tables with >200 runs are processed in batches as 
# specified by the -n flag and $runs_per_batch variable.
# 
# Flags:
# -w : Path to the workflow table.
# -p : Path to the parameter write .R script.
# -g : Path to .R script to gather results into new results table.
# -d : Path to the data directory containing input data objects.
# -r : Path to the rscripts directory.
# -e : Path of the results directory where run outputs are published.
# -t : Regex character string stem of the results table filename.
# -n : Maximum number of runs per batch.
# -c : Whether to clean up files produced by the NextFlow run (T/F).
# 
# Returns:
# Runs the workflow, saving a new results table with unique timestamped filename on
# run success.
#

#------------
# parse flags
#------------
# initialize flag variables
workflow_table_filename=wt.csv
params_write_rscript_filename=r-nf_write-params.R
results_gather_rscript_filename=r-nf_gather-results.R
data_directory=data
rscript_directory=rscript
results_directory=results
results_table=results_table_*
runs_per_batch=200
cleanup=T

#---------------------
# parse provided flags
#---------------------
while getopts ":w:p:g:d:r:e:n:c": name; do
    case "$name" in
        'w')    workflow_table_filename=$OPTARG;;
        'p')    params_write_rscript_filename=$OPTARG;;
		'g')    results_gather_rscript_filename=$OPTARG;;
		'd')	data_directory=$OPTARG;;
		'r')	rscript_directory=$OPTARG;;
		'e')	results_directory=$OPTARG;;
		't')	results_table=$OPTARG;;
		'n')	runs_per_batch=$OPTARG;;
		'c')	cleanup=$OPTARG;;
        \?)     echo "Invalid option provided: -$OPTARG" >&2;;
    esac
done

#-----------------
# manage arguments
#-----------------
workflow_table_path=./$data_directory/$workflow_table_filename
write_param_script_path=./$rscript_directory/$params_write_rscript_filename
results_gather_script_path=./$rscript_directory/$results_gather_rscript_filename

#------------
# check paths
#------------
if ! test -f $workflow_table_path ; then
    echo 'ERROR: workflow table not found at provided path: '$workflow_table_path
    exit 1
fi

if ! test -f "$write_param_script_path"; then
    echo 'ERROR: write parameters .R script not found at provided path: '$write_param_script_path
    exit 1
fi

if ! test -f "$results_gather_script_path"; then
    echo 'ERROR: gather results .R script not found at provided path: '$results_gather_script_path
    exit 1
fi

#----------------------------------
# parse preliminary cleanup options
#----------------------------------
if test -d "$results_directory"; then
    echo 'Removing existing results directory...'
    rm -r $results_directory
fi

echo 'Removing any existing results tables...'
find . -name 'results_table_*' -delete

#-----------------
# check runs count
#-----------------
# get total runs
run_count=$(cat $workflow_table_path | wc -l)
run_count=$(echo "$run_count-1" | bc)

# use generalized algebraic notation to round up batch count
batch_count=$(echo "$run_count+$runs_per_batch-1" | bc)
batch_count=$(echo "$batch_count/$runs_per_batch" | bc)
echo 'Found '$run_count' runs. Parsing in '$batch_count' batches...'

#------------
# handle runs
#------------
batches_remaining=$batch_count
read_start=1
read_end=$(echo "$read_start + $runs_per_batch - 1" | bc)
while (( $batches_remaining > 0 ));
do
	echo $batches_remaining" batches remaining to be processed..."
	echo "updating params.config..."
	echo "read index end is "$read_end
	Rscript $write_param_script_path -f $workflow_table_path -s $read_start -e $read_end
	echo "running workflow..."
	nextflow run main.nf
	echo "Finished batch. Batches left: "$batches_remaining
	batches_remaining=$(echo "$batches_remaining-1" | bc)
	read_start=$(echo "$read_start + $runs_per_batch" | bc)
	read_end=$(echo "$read_start + $runs_per_batch - 1" | bc)
done

#-------------------
# gather new results
#-------------------
echo "gathering results table..."
Rscript $results_gather_script_path

echo "copying results table to results folder..."
rfilename=$(find . -name 'results_table_*')
mv $rfilename ./data/$rfilename

#-------------------
# clean up run files
#-------------------
if [ $cleanup == T ]; then
	echo "Performing cleanup..."
	rm -r 'work'
	rm -r '.nextflow'
	find . -name '*.nextflow.log.*' -delete
	find . -name '.nextflow' -delete
	# find . -name 'results' -delete
else
	echo "Skipping cleanup."
fi

#-------------
# end messages
#-------------
echo "Workflow run success. Returning."