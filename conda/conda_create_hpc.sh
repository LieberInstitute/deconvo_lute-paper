#!/usr/bin bash

# Author: Sean Maden
#
# Make HPC conda environment for manuscript.
#
#
#

srun --pty --x11 --mem=10G bash

module load conda
conda create --name "deconvo_lute-paper" -c bioconda conda-forge r-base=4.3.1
conda activate deconvo_lute-paper

conda install conda-forge::r-biocmanager
conda install bioconda::bioconductor-scuttle=1.12.0
conda install bioconda::bioconductor-scran=1.30.0
conda install conda-forge::dplyr=1.1.3
conda install conda-forge::r-ggplot2=3.4.4
conda install conda-forge::r-gridextra=2.3
conda install conda-forge::ggforce=0.4.1 
conda install bioconda::glmGamPoi=1.14.0 
conda install bioconda::bioconductor-multiassayexperiment=1.28.0
conda install bioconda::bioconductor-summarizedexperiment=1.32.0
conda install bioconda::bioconductor-singlecellexperiment=1.24.0
conda install bioconda::sva=3.50.0
conda install bioconda::bioconductor-limma=3.58.1

conda env export > deconvo_lute-paper.yml