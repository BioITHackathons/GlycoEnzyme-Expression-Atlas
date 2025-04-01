library(dplyr)
library(data.table)

hydrolase <- fread("human_protein_glycohydrolase.csv")
transferase <- fread("human_protein_glycosyltransferase.csv")

glycogene <- fread("human_protein_glycogenes.csv")

setnames(glycogene, "gene_name", "gene_symbol")

hydrolase_genes <- hydrolase$gene_symbol
transferase_genes <- transferase$gene_symbol
glycogene_genes <- glycogene$gene_symbol

all_genes <- c(hydrolase_genes, transferase_genes, glycogene_genes)
unique_genes <- unique(all_genes)

write.table(unique_genes, 
            file = "all_glycogenes.txt", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)