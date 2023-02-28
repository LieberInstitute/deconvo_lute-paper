#!/bin/bash

# Author: Sean Maden
#
# Parse/iterate on a list of workflow tables, copying results tables to data between iterations.
#

#params
datapath=data
rscriptpath=./sh/r-nf_parse-wt-iterations.sh

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
echo "Finished with all workflow table iterations."
