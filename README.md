# IRF6 EMT Analysis

This repository contains code, processed data, and figures used in the analysis of IRF6 as a promoter of an epithelial phenotype.

## Overview
Code and processed data used to generate figures for the manuscript titled “Systems-level analysis identifies IRF6 as an inhibitor of epithelial–mesenchymal transition”.

## Repository structure
- `code/R/`  
  R scripts for correlation analysis and figure generation  
  - `Fig.1.R` generates all panels of Figure 1  
    - Panels A and B use data downloaded via the DepMap (CCLE) portal  
    - Panel C uses data downloaded from the UCSC Xena browser (TCGA)  
    - Panel D uses data downloaded from the respective GEO (GSE) pages  

- `code/python/`  
  Python scripts for generating all panels of Figure 4  
  - `Fig4a.py` can be run directly to generate panel A  
  - `Fig4b.py` requires topology file prefixes (without extensions) as command-line arguments, e.g.:
    ```
    python3 Fig4b.py TS_1 TS_2 TS_3 TS_over_IRF6_1 TS_over_IRF6_2 TS_over_IRF6_3 TS_over_ELF3_1 TS_over_ELF3_2 TS_over_ELF3_3 TS_over_KLF4_1 TS_over_KLF4_2 TS_over_KLF4_3
    ```
  - `Fig4c_d.py` generates panels C and D  
  - `Fig4e_f.py` generates panels E and F  
  - All Figure 4 scripts use data derived from RACIPE simulations  
  - `SuppFig4.py` generates all Supplementary Figure 4 panels using TSV/CSV files in `./data/`

- `data/`  
  Processed TSV/CSV files used for analysis and figure generation

## Data sources
Public datasets including TCGA (via UCSC Xena), CCLE (DepMap), and GEO (NCBI).  
No raw patient-level data is included in this repository.

## Notes
Scripts are provided as used in the analysis and correspond directly to the reported figures.

