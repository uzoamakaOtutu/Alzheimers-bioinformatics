
#-----------------------------------------------
# Volcano plots 
#-----------------------------------------------

# Prepare variables for volcano plots with significant gene results from limma
Volcano_1 <- deg_meta_1
Volcano_2 <- deg_meta_2
Volcano_3 <- deg_meta_3

# Set groups for each dataset
group1 <- factor(sample_info$condition)
group2 <- factor(sample_info_2$condition)
group3 <- factor(sample_info_3$condition)

# Creating a function for plotting volcano plots
volcano_plot <- function(data, group, main_title=NULL, legend_title="Padj < 0.05") {
  # Set colors
  old.pal <- palette(c("#00BFFF", "#FF3030"))
  
  # Set margin size
  par(mar = c(4, 4, 2, 1), cex.main = 0.8)
  
  # Generate title if not specified
  if (is.null(main_title)) {
    main_title <- paste(levels(group)[1], "vs", levels(group)[2])
  }
  print(main_title)
  
  # Plot values
  plot(data$logFC, -log10(data$adj.P.Val), main = main_title,
       xlab = "log2FC", ylab = "-log10(Padj)", pch = 20, cex = 0.4)
  
  # Highlight significant points
  with(subset(data, adj.P.Val < 0.05 & abs(logFC) >= 0.5),
       points(logFC, -log10(adj.P.Val), pch = 20, 
              col = (sign(logFC) + 3) / 2, cex = 0.4))
  
  # Add legend
  legend("bottomleft", title = legend_title,
         legend = c("down", "up"), pch = 20, col = 1:2)
  
  # Restore original palette
  palette(old.pal)
}

# Call function for each dataset to plot volcano plot
volcano_plot(Volcano_1, group1)
volcano_plot(Volcano_2, group2)
volcano_plot(Volcano_3, group3)

# Volcano Plot for GSE5281
png("results/figures/volcano_GSE5281.png", width = 800, height = 600)
volcano_plot(Volcano_1, group1, main_title = "GSE5281")
dev.off()

# Volcano Plot for GSE48350
png("results/figures/volcano_GSE48350.png", width = 800, height = 600)
volcano_plot(Volcano_2, group2, main_title = "GSE48350")
dev.off()

# Volcano Plot for GSE36980
png("results/figures/volcano_GSE36980.png", width = 800, height = 600)
volcano_plot(Volcano_3, group3, main_title = "GSE36980")
dev.off()