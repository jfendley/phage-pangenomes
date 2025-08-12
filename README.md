
This pipeline produces the analysis in J. M. Fendley et. al: "Synteny and linkage decay in bacteriophage pangenomes". 

To recreate the analysis in the paper, first ensure the data is in the correct format. Data from the ActinoBacteriophage Database at https://phagesDB.org can be downloaded using the [phage-download pipeline](https://github.com/jfendley/phage-download) which ensures the correct format. 

Once the data is in the correct format, run the exectuable:
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

The file config/cluster_to_lifestyle.tsv was manually created from https://phagesdb.org/clusters/. The file config/groups_to_process.json was also manually created; see run_everything.sh for more details. There are step-by-step instructions included in run_everything.sh. 

Alternatively, this pipeline is also adaptable for use with other sets of data provided in the same format. Some minor things were coded specifically for the exact dataset used and would need to be modified for a new dataset.

