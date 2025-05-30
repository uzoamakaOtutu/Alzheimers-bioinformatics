
#---------------------------------------------
# Functional enrichment Anlysis using ORA
#---------------------------------------------

# Access gene symbols in a column named "Gene.Symbol"
gene_list <- significant_genes$Gene.Symbol

# Perform ORA - Gene Ontology
ora_results <- enrichGO(
  gene          = gene_list,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL", 
  ont           = "BP",      # "BP" for Biological Process, 
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

# View results
head(ora_results)

# Visualize ORA results
dotplot(ora_results)

# Save GO enrichment results table
go_results_df <- as.data.frame(ora_results)
write.csv(go_results_df, "results/tables/GO_enrichment_results.csv", row.names = FALSE)

# Save GO dotplot
go_dot <- dotplot(ora_results)
ggsave("results/figures/GO_enrichment_dotplot.png", plot = go_dot, width = 8, height = 6)


# Perform ORA - Kegg Pathway
# Convert gene symbols to Entrez IDs
gene_entrez <- bitr(gene_list, fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)

kegg_enrich <- enrichKEGG(gene         = gene_entrez$ENTREZID,
                          organism     = 'hsa', 
                          keyType      = 'kegg', 
                          pvalueCutoff = 0.05)

dotplot(kegg_enrich, showCategory=20)

# Save KEGG enrichment results table
kegg_results_df <- as.data.frame(kegg_enrich)
write.csv(kegg_results_df, "results/tables/KEGG_enrichment_results.csv", row.names = FALSE)

# Save KEGG dotplot
kegg_dot <- dotplot(kegg_enrich, showCategory = 20)
ggsave("results/figures/KEGG_enrichment_dotplot.png", plot = kegg_dot, width = 8, height = 6)

# --------------------------------------------------------
message("GO and KEGG enrichment analysis complete. Results saved.")


#--------------------------------------------------------
# END OF CODE ANALYSIS
#--------------------------------------------------------