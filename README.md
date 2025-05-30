# Alzheimer's Bioinformatics Project

## Overview
Bioinformatics analysis of gene expression, biomarkers, pathway enrichment, and potential therapeutic targets in Alzheimer’s disease using publicly available datasets.

## Features
- Differential gene expression analysis between Alzheimer’s and control samples
- Integration of multiple datasets (GEO)   
- Meta-analysis 
- Gene ontology and pathway enrichment  
- Visualization with volcano plots and enrichment dotplots  


## Prerequisites
- R version 4.0 or higher recommended  
- Internet connection to download datasets (GEO, TCGA)  

## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/uzoamakaOtutu/Alzheimers-bioinformatics.git
cd Alzheimers-bioinformatics
```

### 2. Install Required R Packages
```r
install.packages(c("BiocManager", "ggplot2", "dplyr", "openxlsx"))
BiocManager::install(c("GEOquery", "affy", "limma", "Biobase", "clusterProfiler", "org.Hs.eg.db"))
```

## Running the Pipeline
Open R and run the master script:
```r
source("alz_analysis_master_script.R")
```
This executes all major analysis steps in sequence.

## Project Structure
```
├── data/
│   └── raw/                        # Raw CEL files (excluded via .gitignore)
├── results/
│   ├── figures/                    # Output figures: volcano plots, dotplots
│   └── tables/                     # CSV/XLSX tables with gene lists and enrichment data
├── scripts/                        # Modular analysis scripts
│   ├── data_download_and_preprocessing.R
│   ├── differential_expression.R
│   ├── meta_analysis.R
│   ├── enrichment_analysis.R
│   └── visualizations.R
├── alz_analysis_master_script.R                # Orchestrates full pipeline
└── .gitignore                     # Ignores large raw data files
```

## Notes
- Make sure to run the script in the root directory so relative paths resolve correctly.
- The raw datasets in `data/raw/` are excluded from GitHub with `.gitignore` due to size.

## License
MIT License. See `LICENSE` file for details.
