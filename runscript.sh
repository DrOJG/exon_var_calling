#!/bin/bash -l
# Request 8 hrs of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=8:00:00

# Request 8 gigabyte of RAM per core
#$ -l mem=8G

# Set the name of the job.
#$ -N Exon_align_smk

# Request 12 cores.
#$ -pe smp 12

### Set Working Directory - CHANGE THIS EACH TIME! ###
#$ -wd /home/sejjojg/Scratch/workspace/exon_var_calling/

# Load modules

. $MODULESHOME/init/bash
module load apptainer
module load python/miniconda3/24.3.0-0

source $UCL_CONDA_PATH/etc/profile.d/conda.sh

# Establish snakemake env

conda init bash

conda env create -f ./workflow/envs/snakemake.yml

conda activate snakemake
conda config --set channel_priority strict

# Run
snakemake --sdm conda
