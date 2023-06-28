# Sex shapes cell-type-specific transcriptional signatures of stress exposure in the mouse hypothalamus
## Brivio et al., 2023 Cell Reports
![image](https://github.com/EBrivio/sex_stress/assets/58255155/7290a886-3f4b-4dff-9b70-1c7d52eb9a17)

This repository contains the scripts for the original analysis from the manuscript.

# Web App
This dataset is freely explorable at the following [link](https://male-female-stress.weizmann.ac.il/shinyApp/)

# Content structure
<pre>
Sex x Stress repository  
├── README.md 
└── code
    ├── scRNA-seq analyses (R scripts)  
    │    ├── Whole dataset  
    │    │   ├── Ambient_RNA.R                    # Ambient RNA calculation  
    │    │   ├── Neurons_reclustring.R            # Neurons reclustering  
    │    │   └── Distance_backgrounds.R           # Distance of ARS response between backgrounds   
    │    └── Oligodendrocytes analysis  
    │         ├── OL_DEGs_clustering.Rmd          # Genes clustering  
    │         ├── OL_pseudotime.R                 # Pseudotime analysis 
    │         └── OL_CCInx.R                      # CCInx permutation analysis 
    │  
    └── RNAscope analysis (python and groovy scripts)  
        ├── RNAscope_processing_part1_images_preparation.groovy         # Images preparation (groovy)  
        ├── RNAscope_processing_part2_images_preprocessing.ipynb        # Images preprocessing (py)
        ├── RNAscope_processing_part3_roi_expansion_nooverlap.groovy    # ROI espansion (groovy)  
        └── RNAscope_processing_part4_signal_identification.ipynb       # Signal detection (py)
</pre>

