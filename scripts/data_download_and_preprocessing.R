
#-----------------------------------------------------
# Download and process raw data for each data set
#-----------------------------------------------------

#--------------------
# DATA 1 - GSE5281
#--------------------
# Downloading supplementary files for the GEO data set
#getGEOSuppFiles("GSE5281")
#untar("GSE5281_RAW.tar", exdir = 'data5281/')

# Reading the raw CEL files into an AffyBatch object
cel_files_path_1 <- "data/raw/data5281"
cel_files_1 <- list.files(cel_files_path_1, "\\.CEL\\.gz$", full.names = TRUE)
lapply(cel_files_1, gunzip, remove = TRUE, overwrite = TRUE)
raw_data1 <- ReadAffy(celfile.path = cel_files_path_1)

# Normalizing the raw data using the RMA method
normalized_data1 <- rma(raw_data1)

# Extracting expression values from normalized data
exprs_data_1 <- exprs(normalized_data1)

#--------------------
# DATA 2 - GSE48350
#--------------------

# Downloading supplementary files for the GEO data set
#getGEOSuppFiles("GSE48350")
#untar("GSE48350_RAW.tar", exdir = 'data48350/')

# Reading the raw CEL files into an AffyBatch object
cel_files_path_2 <- "data/raw/data48350"
cel_files_2 <- list.files(cel_files_path_2, "\\.CEL\\.gz$", full.names = TRUE)
lapply(cel_files_2, gunzip, remove = TRUE, overwrite = TRUE)
raw_data2 <- ReadAffy(celfile.path = cel_files_path_2)

# Normalizing the raw data using the RMA method
normalized_data2 <- rma(raw_data2)

# Extracting expression values from normalized data
exprs_data_2 <- exprs(normalized_data2)

#--------------------
# DATA 3 - GSE36980
#--------------------

# Downloading supplementary files for the GEO data set
#getGEOSuppFiles("GSE36980")
#untar("GSE36980_RAW.tar", exdir = 'data36980/')

# Define the directory where CEL files are located
cel_files_path <- "data/raw/data36980/"
library(pd.hugene.1.0.st.v1)
library(oligo)
# Create a list of CEL files in the specified directory
cel_files <- list.files(cel_files_path, 
                        pattern = "\\.CEL\\.gz$", 
                        full.names = TRUE)

# Read the CEL files into an oligo object
raw_data3 <- read.celfiles(cel_files)

# Normalizing the raw data using the RMA method
normalized_data3 <- rma(raw_data3)

# Extracting expression values from normalized data
exprs_data_3 <- exprs(normalized_data3)

#----------------------------------------------
# EXTRACTING SAMPLE DATA FOR CLEANING
#----------------------------------------------

#Data 1 - GSE5281
# Extracting sample information using GEOquery
gse_1 <- getGEO("GSE5281", GSEMatrix = TRUE)
info_data_1 <- gse_1[[1]]
sample_data_1 <- pData(info_data_1)

#Data 2 - GSE48350
# Extracting sample information using GEOquery
gse_2 <- getGEO("GSE48350", GSEMatrix = TRUE)
info_data_2 <- gse_2[[1]]
sample_data_2 <- pData(info_data_2)

#Data 3 - GSE36980
# Extracting sample information using GEOquery
gse_3 <- getGEO("GSE36980", GSEMatrix = TRUE)
info_data_3 <- gse_3[[1]]
sample_data_3 <- pData(info_data_3)

#-----------------------------------------------------
# SAMPLE DATA CLEANING FOR DATA 1
#-----------------------------------------------------

# Displaying the first few rows to understand its structure
head(sample_data_1)

# Renaming columns
colnames(sample_data_1)[colnames(sample_data_1) == 
                          "characteristics_ch1.8"] <- "condition"

# Keeping columns of interest
columns_to_keep <- c("condition")

# Creating a new data frame with only the specified columns
sample_data_clean <- sample_data_1[, columns_to_keep, drop = FALSE]

# Removing prefixes "Disease State: " and "Sex: "
sample_data_clean$condition <- gsub("(?i)^disease state: ", "", 
                                    sample_data_clean$condition, perl = TRUE)

# Modifying the 'condition' columns
sample_data_clean$condition <- gsub("normal", "control", 
                                    sample_data_clean$condition, ignore.case = TRUE)

sample_data_clean$condition <- gsub("Alzheimer's Disease", "AD", 
                                    sample_data_clean$condition, ignore.case = TRUE)

# Replacing non-breaking spaces (U+00A0) with regular spaces
sample_data_clean$condition <- gsub("\u00A0", " ", 
                                    sample_data_clean$condition)

# Trimming any leading and trailing whitespace characters
sample_data_clean$condition <- trimws(sample_data_clean$condition)

# Assigning cleaned data to sample_info
sample_info <- sample_data_clean

#-------------------------------------
# SAMPLE DATA CLEANING FOR DATA 2
#-------------------------------------

# Displaying the first few rows to understand its structure
head(sample_data_2)

# Renaming columns
colnames(sample_data_2)[colnames(sample_data_2) == 
                          "source_name_ch1"] <- "condition"

# Keeping columns of interest
columns_to_keep_2 <- c("condition")

# Creating a new data frame with only the specified columns
sample_data_clean_2 <- sample_data_2[, columns_to_keep_2, drop = FALSE]

# Replacing with "control" or "AD" based on "_AD" presence
sample_data_clean_2$condition <- ifelse(
  grepl("_AD", sample_data_clean_2$condition), "AD2", "control"
)

# Assigning cleaned data to sample_info
sample_info_2 <- sample_data_clean_2

#-------------------------------------
# SAMPLE DATA CLEANING FOR DATA 3
#-------------------------------------

# Displaying the first few rows to understand its structure
head(sample_data_3)

# Renaming columns
colnames(sample_data_3)[colnames(sample_data_3) == "title"] <- "condition"

# Keeping columns of interest
columns_to_keep_3 <- c("condition")

# Creating a new data frame with only the specified columns
sample_data_clean_3 <- sample_data_3[, columns_to_keep_3, drop = FALSE]

# Replacing with "control" or "AD" based on "_AD" presence
sample_data_clean_3$condition <- ifelse(
  grepl("non-AD", sample_data_clean_3$condition), "control", "AD3"
)

# Assigning cleaned data to sample_info
sample_info_3 <- sample_data_clean_3
