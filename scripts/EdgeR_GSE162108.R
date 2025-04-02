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

write.csv(res_glyco, file = "GSE162108_DE_results.csv", row.names = TRUE)