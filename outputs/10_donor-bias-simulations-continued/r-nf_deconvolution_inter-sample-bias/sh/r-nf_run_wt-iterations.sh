#!/bin/bash

# Author: Sean Maden
#
# Parse/iterate on a list of workflow tables, copying results tables to data between iterations.
#

#params
datapath=data
rscriptpath=./sh/r-nf_parse-wt-iterations.sh
rscript_directory=rscript
results_gather_rscript_filename=r-nf_gather-results.R
results_gather_script_path=./$rscript_directory/$results_gather_rscript_filename

# detect tables
workflow_table_filename=`find $datapath -type f | grep 'workflow-table-iter*'`
wt=($workflow_table_filename)
wtlen=`echo "${#wt[@]}"`
echo "Found "$wtlen" workflow tables for iterations..."

echo "Doing iterations..."
for filename in $workflow_table_filename
do
	echo "Working on workflow table "$filename"..." 
	basename="$(basename -- $filename)"
	bash $rscriptpath -w $basename
	echo "Run complete. Continuing..."
done
echo "Workflow run success."

#-------------------
# gather new results
#-------------------
echo "gathering results table..."
Rscript $results_gather_script_path

echo "copying results table to results folder..."
rfilename=$(find . -name 'results_table_*')
cp $rfilename $rfilename

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
echo "Finished with all workflow table iterations."
