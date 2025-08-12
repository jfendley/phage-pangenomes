# Phage Pangenome Analysis Pipeline

This pipeline produces the analysis in
> J. M. Fendley et. al\
> _"Synteny and linkage decay in bacteriophage pangenomes"._

The pipeline is designed to process bacteriophage genomic data to reproduce all figures from the study. It can also be adapted for other datasets.

## Running the Pipeline

The pipeline analyzes data downloaded from the ActinoBacteriophage Database (https://phagesDB.org). As a preliminary step, you can use the [phage-download pipeline](https://github.com/jfendley/phage-download) to download and prepare the data in the correct format.

Once the data has been downloaded, run the executable:
```
bash run_everything.sh
```
Configured for 64GB memory and 20 cores this will take approximately 14 hours.

You can consult the `run_everything.sh` file for more instructions and details on how to run the pipeline.

The files `config/cluster_to_lifestyle.tsv` and `config/groups_to_process.json` were manually created from https://phagesdb.org/clusters/, and are adapted for the data used in the study. They can be modified in case you want to run the pipeline on a different set of data, see `run_everything.sh` for more details.

## Results

The output will be created in the `results` directory. The directory has the following structure, with all of the figures in the paper in the folder `paper_figures`.
```
└── results/
    ├── groups/
    │   ├── E/
    │   │   ├── E_core/
    |   |   ├── E_accessory/
    │   │   ├── E_defense/
    │   │   └── E_phams/
    │   └── etc...
    ├── all_groups_figures/
    ├── all_groups_files/
    ├── defense_finder/
    └── paper_figures/
```
