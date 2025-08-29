# Phage Pangenome Analysis Pipeline

This pipeline produces the analysis in
> J. M. Fendley, M. Molari, R. A. Neher, and B. I. Shraiman\
> _"[Synteny and linkage decay in bacteriophage pangenomes](https://www.biorxiv.org/content/10.1101/2025.08.12.669904v1)"._

The pipeline is designed to process bacteriophage genomic data to reproduce all figures from the study. It can also be adapted for other datasets.

## Running the Pipeline

The pipeline analyzes data downloaded from the Actinobacteriophage Database (https://phagesDB.org). As a preliminary step, you can use the [phage-download pipeline](https://github.com/jfendley/phage-download) to download and prepare new data in the correct format.

To use the exact data used in the paper, extract `data.zip`:
```
unzip data.zip
```

Once the data has been extracted or new data has been downloaded, run the executable:
```
bash run_everything.sh
```
Configured for 64GB memory and 20 cores this will take approximately 14 hours with the data used in the paper.

You can consult the `run_everything.sh` file for more instructions and details on how to run the pipeline.

The file `config/cluster_to_lifestyle.tsv` was manually created from https://phagesdb.org/clusters/, and the file `config/groups_to_process.json` was manually created from inspection of results. They, and some clearly specified parts of scripts, should be modified if the pipeline is run on a different set of data; see `run_everything.sh` for more details.

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
