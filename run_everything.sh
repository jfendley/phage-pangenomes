#!/bin/bash


# the following has been tested on 64GB RAM computers. Some of the 
#   analyses (specifically the linkage disequilibrium analyses) require
#   large amounts of memory. The code could be modified (e.g. delete no 
#   longer used variables, calculate linkage only at short distances) to 
#   reduce memory usage if needed.

# first, ensure the data is is in the same format as the output of the phage-
#   download pipeline

# create and activate the conda environment. creating the conda environment 
#   and updating the DefenseFinder models requires internet

conda env create -f environment.yml
source activate base 
conda activate phage-pangenomes

# update defense-finder
defense-finder update --models-dir defense_finder_models

# now the rest of the script does not require internet access. However, LaTeX
#   and TeX Live need to be installed to recreate the paper figures. If they 
#   are not installed, you can modify parameters.py as suggested in the script 
#   for it to run without LaTeX.

# run the pipeline that checks core synteny, and creates sytenic groups 
#   the number of cores can be modified based on computing resources
#   with 20 cores (64GB memory) this step takes approximately 15 minutes

snakemake --cores 20 -s snakefile_synteny all

# then run all of the non-linkage analyses (about 35 minutes)

snakemake --cores 20 -s snakefile_analysis all_except_linkage

# then run the linkage analyses. this uses a lot of memory so is best on
#   only one core, and takes about nine and a half hours

snakemake --cores 1 -s snakefile_analysis linkage

# then the file config/groups_to_process.json was manually created by 
#   inspecting the results. If new data is used, a new file would be needed.

# then run the non-linkage subgroup + processing analyses
#   The following step should only take a minute or so. 

snakemake --cores 20 -s snakefile_subgroups all_except_linkage

# finally, run the subgroup + processing linkage analyses
#   this step takes around three hours. 

snakemake --cores 1 -s snakefile_subgroups linkage
