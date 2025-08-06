
This pipeline reproduces the analysis in J. M. Fendley et. al: "Synteny and linkage decay in bacteriophage pangenomes". 

To recreate the figures and all of the analysis using the data provided, download the data.zip file, add it to this folder, and then run the exectuable:
```
bash run_everything.sh
```
Configured for 64GB memory and 20 cores this will take approximately 14 hours. 

The results folder will be in the following format, with all of the figures in the paper in the folder paper_figures. Figure 2 and parts of Figure 3 and 4 were manually drawn using Inkscape.
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

Alternatively, this pipeline is also adaptable for use with other sets of data provided in the same format. New data from the AcinoBacteriophage Database at https://phagesDB.org can be downloaded using the [phage-download pipeline](https://github.com/jfendley/phage-download)
 which ensures the correct format. 

The file config/cluster_to_lifestyle.tsv was manually created from https://phagesdb.org/clusters/. The file config/groups_to_process.json was also manually created; see run_everything.sh for more details. There are step-by-step instructions included in run_everything.sh.

