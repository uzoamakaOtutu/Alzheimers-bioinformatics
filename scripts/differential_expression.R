
#-------------------------------------------------------------
# Differential expression analysis using Limma
#-------------------------------------------------------------
#---------------------------------------
# DATA 1 - Alzheimer's vs. Control - (GSE5281)
#---------------------------------------

# Creating a group factor based on the 'condition' column
group <- factor(sample_info$condition)
levels(group)

# Designing a matrix for the model
design <- model.matrix(~ 0 + group)

# Checking the column names in the design matrix
print(colnames(design))

# Fitting the linear model to the normalized data
fit <- lmFit(exprs_data_1, design)

# Creating a contrast matrix for comparisons (AD vs control)
contrast.matrix <- makeContrasts(
  AD_vs_Control = groupAD - groupcontrol, 
  levels = design
)

# Applying the contrasts to the fit object
fit2 <- contrasts.fit(fit, contrast.matrix)

# Computing empirical Bayes statistics
fit2 <- eBayes(fit2)

# Extracting the top differentially expressed genes
results <- topTable(fit2, adjust = "BH", number = Inf)

# View the results
print(results)

#----------------------------------------------
# DATA 2 - Alzheimer's vs. Control (GSE48350)
#----------------------------------------------

# Creating a group factor based on the 'condition' column
group <- factor(sample_info_2$condition)
levels(group)
print(group)

# Designing a matrix for the model
design <- model.matrix(~ 0 + group)

# Checking the column names in the design matrix
print(colnames(design))

# Fitting the linear model to the normalized data
fit <- lmFit(exprs_data_2, design)

# Creating a contrast matrix for comparisons (AD2 vs control)
contrast.matrix <- makeContrasts(
  AD2_vs_Control = groupAD2 - groupcontrol, 
  levels = design
)

# Applying the contrasts to the fit object
fit2 <- contrasts.fit(fit, contrast.matrix)

# Computing empirical Bayes statistics
fit2 <- eBayes(fit2)

# Extracting the top differentially expressed genes
results_2 <- topTable(fit2, adjust = "BH", number = Inf)

# View the results
print(results_2)

#----------------------------------------------
# DATA 3 - Alzheimer's vs. Control (GSE36980)
#----------------------------------------------

# Creating a group factor based on the 'condition' column
group <- factor(sample_info_3$condition)
levels(group)
print(group)

# Designing a matrix for the model
design <- model.matrix(~ 0 + group)

# Checking the column names in the design matrix
print(colnames(design))

# Fitting the linear model to the normalized data
fit <- lmFit(exprs_data_3, design)

# Creating a contrast matrix for comparisons (AD3 vs control)
contrast.matrix <- makeContrasts(
  AD3_vs_Control = groupAD3 - groupcontrol, 
  levels = design
)

# Applying the contrasts to the fit object
fit2 <- contrasts.fit(fit, contrast.matrix)

# Computing empirical Bayes statistics
fit2 <- eBayes(fit2)

# Extracting the top differentially expressed genes
results_3 <- topTable(fit2, adjust = "BH", number = Inf)

# View the results
print(results_3)