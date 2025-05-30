
#----------------------------------------------------------------
# Extracting Significant genes - Upregulation and Downregulation
#----------------------------------------------------------------

# Function to filter significant genes 
# and categorize into upregulated and downregulated
get_sig_genes <- function(deg_results, adj_pval_threshold = 0.05, 
                          logFC_threshold = 0.5) {
  # Filter for significant genes based on adjusted p-value 
  #and absolute log fold change
  significant_genes <- deg_results[
    deg_results$adj.P.Val < adj_pval_threshold & 
      abs(deg_results$logFC) > logFC_threshold, 
  ]
  
  # Order by adjusted p-value
  significant_genes <- significant_genes[
    order(significant_genes$Gene.Symbol, significant_genes$adj.P.Val), 
  ]
  
  # For each gene, keep the entry with the lowest adjusted p-value
  unique_genes <- significant_genes[!duplicated(significant_genes$Gene.Symbol), ]
  
  # Sort by adjusted p-value
  sorted_genes <- unique_genes[order(unique_genes$adj.P.Val), ]
  
  return(sorted_genes)
}

# Call the function with default thresholds
#Getting significant genes for Data 1 - AD vs Control
gene_dfs <- get_sig_genes(deg_data_f)

#Getting significant genes for Data 2 - AD vs Control
gene_dfs2 <- get_sig_genes(deg_data_f_2)

#Getting significant genes for Data 3 - AD vs Control
gene_dfs3 <- get_sig_genes(deg_data_f_3)


# Get the top 50 genes for each dataset
top_genes_df1 <- gene_dfs[1:50, ]
top_genes_df2 <- gene_dfs2[1:50, ]
top_genes_df3 <- gene_dfs3[1:50, ]


# Define logFC threshold
logfc_threshold <- 0.5

# Filter upregulated and downregulated genes based on logFC
# For Dataset 1
upreg_genes_df1 <- gene_dfs[gene_dfs$logFC > logfc_threshold, ]
downreg_genes_df1 <- gene_dfs[gene_dfs$logFC < -logfc_threshold, ]

# For Dataset 2
upreg_genes_df2 <- gene_dfs2[gene_dfs2$logFC > logfc_threshold, ]
downreg_genes_df2 <- gene_dfs2[gene_dfs2$logFC < -logfc_threshold, ]

# For Dataset 3
upreg_genes_df3 <- gene_dfs3[gene_dfs3$logFC > logfc_threshold, ]
downreg_genes_df3 <- gene_dfs3[gene_dfs3$logFC < -logfc_threshold, ]

#send data to excel file
wb <- createWorkbook()

# Add each data frame to a separate sheet

# Data frame 1
addWorksheet(wb, "GSE5281_Significant_Genes")
writeData(wb, sheet = "GSE5281_Significant_Genes", gene_dfs)

addWorksheet(wb, "GSE48350_Significant_Genes")
writeData(wb, sheet = "GSE48350_Significant_Genes", gene_dfs2)

addWorksheet(wb, "GSE36980_Significant_Genes")
writeData(wb, sheet = "GSE36980_Significant_Genes", gene_dfs3)

# Save the workbook
saveWorkbook(wb, file = "results/tables/significant_genes_combined.xlsx", overwrite = TRUE)


# Create a dataframe for the top 50 genes of each dataset
combined_top_genes <- data.frame(GSE5281 = top_genes_df1$Gene.Symbol,
                                 GSE48350 = top_genes_df2$Gene.Symbol,
                                 GSE36980 = top_genes_df3$Gene.Symbol)

# Writing the results to a CSV file
write.csv(
  combined_top_genes, 
  file = "results/tables/combined_top_genes_all_dataset.csv", 
  row.names = TRUE
)
