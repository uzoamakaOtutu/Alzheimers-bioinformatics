
#-------------------------------------------
# MAPPING PROBES TO ENTREZ_ID SYMBOLS
#-------------------------------------------
# Define the path to the downloaded gzipped and destination file
#For GPL6244 annotation (GSE36980)
#gz_file <- "data/raw/GPL6244.annot.gz"
txt_file <- "data/raw/GPL6244.annot.txt"

# Unzip the annotation file
#gunzip(gz_file, destname = txt_file, overwrite = TRUE)

# Reading annotation file 2 into environment (For Affy 1.0).
annot_2 <- read.table(txt_file, sep = "\t", header = TRUE, 
                      stringsAsFactors = FALSE, skip = 27, 
                      fill = TRUE, quote = "")
colnames(annot_2)[colnames(annot_2) == "Gene.symbol"] <- "Gene.Symbol"

# Reading annot file into environment(Affy 2.0)
annot <- read.table("data/raw/GPL570-55999.txt", sep = "\t", header = TRUE, 
                    stringsAsFactors = FALSE, skip = 16, 
                    fill = TRUE, quote = "")

# Viewing the first few rows to understand the structure
head(annot)
probe_to_gene <- annot[, c("ID", "Gene.Symbol")]
probe_to_gene2 <- annot_2[, c("ID", "Gene.Symbol", "Platform_SPOTID")]

# Map for AD disease - condition and control - data1 (GSE5281)
results$ProbeID <- rownames(results)
deg_data <- merge(results, probe_to_gene, by.x = "ProbeID", by.y = "ID")

# Map for AD disease condition and control - data2 (GSE48350)
results_2$ProbeID <- rownames(results_2)
deg_data_2 <- merge(results_2, probe_to_gene, by.x = "ProbeID", by.y = "ID")

# Map for AD disease condition and control - data3 (GSE36980)
results_3$ProbeID <- rownames(results_3)
deg_data_3 <- merge(results_3, probe_to_gene2, by.x = "ProbeID", by.y = "ID")

# Excluding control probes for Data 1
deg_data_f <- deg_data[!grepl("^AFFX-", deg_data$ProbeID), ]

# Excluding control probes for Data 2
deg_data_f_2 <- deg_data_2[!grepl("^AFFX-", deg_data_2$ProbeID), ]

# Excluding control probes for Data 3
deg_data_f_3 <- deg_data_3[
  !grepl("control", deg_data_3$Platform_SPOTID, ignore.case = TRUE), 
]

# Remove the Platform_SPOTID column
deg_data_f_3 <- deg_data_f_3[, !names(deg_data_f_3) %in% "Platform_SPOTID"]

# Filtering out rows where Gene.Symbol is empty or NA - Data 1
deg_data_f <- deg_data_f[
  !is.na(deg_data_f$Gene.Symbol) & deg_data_f$Gene.Symbol != "", 
]

# Filtering out rows where Gene.Symbol is empty or NA - Data 2
deg_data_f_2 <- deg_data_f_2[
  !is.na(deg_data_f_2$Gene.Symbol) & deg_data_f_2$Gene.Symbol != "", 
]

# Filtering out rows where Gene.Symbol is empty or NA - Data 3
deg_data_f_3 <- deg_data_f_3[
  !is.na(deg_data_f_3$Gene.Symbol) & deg_data_f_3$Gene.Symbol != "", 
]

remove_duplicates <- function(df) {
  # Order by Gene.Symbol and then by adjusted p-value
  ordered_genes <- df[order(df$Gene.Symbol, df$adj.P.Val), ]
  
  # For each gene, keep the entry with the lowest adjusted p-value
  unique_genes <- ordered_genes[!duplicated(ordered_genes$Gene.Symbol), ]
  
  # Sort by adjusted p-value
  sorted_genes <- unique_genes[order(unique_genes$adj.P.Val), ]
  
  return(sorted_genes)
}

# Apply the function to each dataset
clean_limma_1 <- remove_duplicates(deg_data_f)
clean_limma_2 <- remove_duplicates(deg_data_f_2)
clean_limma_3 <- remove_duplicates(deg_data_f_3)

# Specifying the output file path
#gene_file_1 <- "GSE_5281_gene_results_disease.csv"

#gene_file_2 <- "GSE_48350_gene_results_disease.csv"

#gene_file_3 <- "GSE_36980_gene_results_disease.csv"

# Writing the results to a CSV file
#write.csv(clean_limma_1, file = gene_file_1, row.names = TRUE)

#write.csv(clean_limma_2, file = gene_file_2, row.names = TRUE)

#write.csv(clean_limma_3, file = gene_file_3, row.names = TRUE)

