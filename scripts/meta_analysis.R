
#---------------------------------
#      META ANALYSIS
#---------------------------------
#---------------------------------
# Using FISHER'S METHOD
#---------------------------------
library(metafor)
# Function to rank and filter dataset
rank_and_filter <- function(df) {
  # Order by p-value
  ordered_genes <- df %>%
    arrange(P.Value)
  
  # Remove duplicates by keeping the entry with the lowest p-value
  unique_genes <- ordered_genes[!duplicated(ordered_genes$Gene.Symbol), ]
  
  return(unique_genes)
}

# Apply the function to each dataset
deg_meta_1 <- rank_and_filter(deg_data_f)
deg_meta_2 <- rank_and_filter(deg_data_f_2)
deg_meta_3 <- rank_and_filter(deg_data_f_3)


# Each dataframe should have columns: ProbeID, GeneSymbol, logFC, and p.value

# Merge the data based on common GeneSymbol - disease
deg_combined <- Reduce(function(x, y) merge(x, y, by = "Gene.Symbol", all = TRUE),
                       list(deg_meta_1, deg_meta_2, deg_meta_3))

# Identifying p-value columns
pval_columns <- grep("P.Value", colnames(deg_combined), value = TRUE)

# Ensure the p-values are numeric
deg_combined[pval_columns] <- lapply(
  deg_combined[pval_columns], 
  function(col) as.numeric(as.character(col))
)

# Calculate Fisher's statistic for Disease
deg_combined$fisher_stat <- apply(
  deg_combined[, pval_columns], 1, 
  function(pvals) {
    pvals <- pvals[!is.na(pvals)]  # Remove NAs
    if (length(pvals) == 0) return(NA)  # Handle case with no valid p-values
    -2 * sum(log(pvals))
  }
)

# Check for non-numeric values in fisher_stat and handle them - disease
non_numeric_fisher_stat <- !is.numeric(deg_combined$fisher_stat) | 
  is.na(deg_combined$fisher_stat)

# Check if there are any non-numeric or NA values
if (any(non_numeric_fisher_stat, na.rm = TRUE)) {
  # Output indices of problematic rows
  warning("Non-numeric values detected in fisher_stat at indices: ",
          paste(which(non_numeric_fisher_stat), collapse = ", "))
  # Set problematic rows to NA
  deg_combined$fisher_stat[non_numeric_fisher_stat] <- NA
}

# Remove rows with NA in any p-value columns
pval_columns <- grep("P.Value", colnames(deg_combined), value = TRUE)
deg_combined_cleaned <- deg_combined[
  complete.cases(deg_combined[pval_columns]), 
]

# Calculate degrees of freedom for each row - disease
deg_combined_cleaned$df <- apply(
  deg_combined_cleaned[, pval_columns], 1, 
  function(pvals) {
    2 * sum(!is.na(pvals))  # 2 times the number of non-NA p-values
  }
)

# Ensure df is numeric - disease
deg_combined_cleaned$df <- as.numeric(deg_combined_cleaned$df)

# Calculate p-values from Fisher's statistic - disease
deg_combined_cleaned$fisher_p <- pchisq(
  deg_combined_cleaned$fisher_stat, 
  df = deg_combined_cleaned$df, 
  lower.tail = FALSE
)

# Adjusting p-values for multiple testing - disease
deg_combined_cleaned$adj_fisher_p <- p.adjust(
  deg_combined_cleaned$fisher_p, 
  method = "BH"
)

# Filtering significant genes with adjusted Fisher p-value, 
#absolute logFC >= 0.5, and at least one dataset with average expression > 3
significant_genes <- deg_combined_cleaned %>%
  filter(adj_fisher_p < 0.05 & 
           (abs(logFC) >= 0.5 | abs(logFC.x) >= 0.5 | abs(logFC.y) >= 0.5) & 
           (AveExpr > 3 | AveExpr.x > 3 | AveExpr.y > 3))

# Sort the significant results by adjusted Fisher p-value
sorted_significant_results <- significant_genes[
  order(significant_genes$adj_fisher_p), 
]

# Select the top 50 genes
top_50_genes_fisher <- head(sorted_significant_results, 50)

# Display the top 50 genes
print(top_50_genes_fisher$Gene.Symbol)

#specify file path
#Fisher_meta_results <- "Fisher_meta_results.csv"
#write.csv(
#  sorted_significant_results, 
#  file = Fisher_meta_results, 
#  row.names = TRUE
#)

#specify file path
Fisher_meta_results_top_50 <- "results/tables/Fisher_meta_results_top_50.csv"
write.csv(top_50_genes_fisher, Fisher_meta_results_top_50, row.names = TRUE)
