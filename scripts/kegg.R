kegg1 <- fread("GSE140261_DE_results.csv")
kegg2 <- fread("GSE162108_DE_results.csv")

library(clusterProfiler)
library(org.Hs.eg.db)     
library(enrichplot)
library(DOSE)

gene_symbols <- kegg1$gene_symbol

# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

entrez_ids <- unique(entrez_ids$ENTREZID)

kegg_enrich <- enrichKEGG(gene = entrez_ids,
                          organism = 'hsa',
                          pvalueCutoff = 0.05)

kegg_df <- as.data.frame(kegg_enrich)