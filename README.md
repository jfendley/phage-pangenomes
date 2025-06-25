
This pipeline reproduces the analysis in J. M. Fendley et. al: "Synteny and linkage decay in bacteriophage pangenomes". 

To simply recreate the figures and all of the analysis using the data provided used simply run the exectuable:
```
./run_everything.sh
```
Configured for 64GB memory and 20 cores, this will take approximately 15 hours. 

The results folder is in the following format:
```
├── results/
│   ├── groups/
│   │   ├── E/
│   │   │   ├── E_core/
|   |   |   ├── E_accessory/
│   │   │   ├── E_defense/
│   │   │   └── E_phams/
│   │   └── etc...
│   ├── all_groups_figures/
│   ├── all_groups_files/
│   ├── defense_finder/
│   └── paper_figures/
```

This pipeline is also adaptable to other sets of data provided in the same format. New data from the AcinoBacteriophage
Database at https://phagesDB.org can be downloaded using the phage-download pipeline which ensures the correct format. The original 
data can be accessed by unzipping the data file. 
```
unzip data.zip
```
The file config/cluster_to_lifestyle.tsv was manually created from https://phagesdb.org/clusters/

Once the data is sufficiently formated, create and activate the conda environment.
```
conda env create -f environment.yml
conda activate phage-pangenomes
```

First, run the pipeline that checks core synteny, and creates sytenic groups. The number of cores can be modified based on computing resource. With 20 cores (64GB memory) this step takes approximately 25 minutes.

```
snakemake --cores 20 -s snakefile_synteny all
```

The next step takes about 45 minutes. 
```
snakemake --cores 20 -s snakefile_analysis all_except_linkage
```
This is the time-consuming step, and it takes about 8 hours or so. 
```
snakemake --cores 1 -s snakefile_analysis linkage
```

Then the file config/groups_to_process.json was manually created by inspecting the results. If new data is used, a new file would be needed.

The following step should only take a minute or so. 

```
snakemake --cores 20 -s snakefile_subgroups all_except_linkage
```
And the final step takes
```
snakemake --cores 1 -s snakefile_subgroups linkage
```