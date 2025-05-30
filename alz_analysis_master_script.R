
# -----------------------------------------------------------------------
# Master Script: master_script.R
# Project: Alzheimer's Bioinformatics Analysis
# Author: Uzoamaka Otutu
# Description: Runs the full pipeline by sourcing scripts 
#.             from the 'scripts/' folder in sequence.
# -----------------------------------------------------------------------


# Install packages and load libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GEOquery", "affy", "limma","clusterProfiler","org.Hs.eg.db","Biobase", "metafor", "impute"))

library(GEOquery)
library(affy)
library(R.utils)
library(limma)
library(Biobase)
library(dplyr)
library(openxlsx)
library(impute)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)


#--------------------------------
# Run scripts in order
#--------------------------------

# Step 1: Download and preprocess raw data
message("Step 1: Downloading and preprocessing data...")
source("scripts/data_download_and_preprocessing.R")

# Step 2: Perform differential expression analysis
message("Step 2: Performing differential expression analysis...")
source("scripts/differential_expression.R")

# Step 3: Map gene symbols or Ensembl IDs
message("Step 3: Mapping gene symbols...")
source("scripts/mapping_gene_symbols.R")

# Step 4: Filter significant DEGs (e.g., logFC, p-value)
message("Step 4: Filtering significant DEGs...")
source("scripts/filter_significant_genes.R")

# Step 5: Perform meta-analysis
message("Step 5: Performing meta-analysis...")
source("scripts/meta_analysis.R")

# Step 6: Run enrichment analysis (GO, KEGG)
message("Step 6: Running enrichment analysis...")
source("scripts/enrichment_analysis.R")

# Step 7: Generate visualizations
message("Step 7: Generating visualizations...")
source("scripts/visualizations.R")

# Done
message("All steps completed successfully.")


