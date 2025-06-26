#!/bin/bash

# unzip the data file (alternatively could download new data)

unzip data.zip 

# create and activate the conda environment

conda env create -f environment.yml
source activate base 
conda activate phage-pangenomes

# run the pipeline that checks core synteny, and creates sytenic groups 
#   the number of cores can be modified based on computing resources
#   with 20 cores (64GB memory) this step takes approximately 25 minutes

snakemake --cores 20 -s snakefile_synteny all

# then run all of the non-linkage analyses (about 45 minutes)

snakemake --cores 20 -s snakefile_analysis all_except_linkage

# then run the linkage analyses. this uses a lot of memory so is best on
#   only one core, and takes about eight hours

snakemake --cores 1 -s snakefile_analysis linkage

# then run the non-linkage subgroup + processing analyses
#   The following step should only take a minute or so. 

snakemake --cores 20 -s snakefile_subgroups all_except_linkage

# finally, run the subgroup + processing linkage analyses
#   this step takes around six hours

snakemake --cores 1 -s snakefile_subgroups linkage
