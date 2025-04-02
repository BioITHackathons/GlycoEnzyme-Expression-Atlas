expr <- read.delim("data/GSE162108_Falzarano_raw_counts.txt", row.names = 1, check.names = FALSE)
expr_filtered <- expr[, c("IG-m", "C-m")]

#EDGE R
group <- factor(c("DMD", "Control"))
dge <- DGEList(counts = expr_filtered, group = group)

dge$common.dispersion <- 0.1

et <- exactTest(dge, dispersion = 0.1)

res <- topTags(et, n = Inf)$table

# Filter by the glycogenes
glyco_genes <- readLines("data/all_glycogenes.txt") 

res_glyco <- res[rownames(res) %in% glyco_genes, ]

res_glyco$regulation <- ifelse(res_glyco$logFC > 0, "Up", "Down")

library(clusterProfiler)
library(org.Hs.eg.db)     
library(enrichplot)
library(DOSE)

gene_symbols <- rownames(res_glyco)

# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

entrez_ids <- unique(entrez_ids$ENTREZID)

kegg_enrich <- enrichKEGG(gene = entrez_ids,
                          organism = 'hsa',
                          pvalueCutoff = 0.05)




