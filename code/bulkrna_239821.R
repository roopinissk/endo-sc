BiocManager::install("Rsubread")
library(dplyr)
library(tidyverse)
BiocManager::install("GEOquery")
library(GEOquery)
library(readxl)
library(tidyr)
sessionInfo("dplyr")
library(magrittr)
"
This is the bulk rna dataset
to do with this dataset: 
1. DEGs for immune related genes; 
2. gene ontology enrichment analysis for immune  functions
3. pathway analysis to identify immune signaling pathways activated in blastocysts
"

input_dir <- "~/Documents/nyu/spring/ms_thesis/GSE_239821/gse239821_data"

output_dir <- "~/Documents/nyu/spring/ms_thesis/GSE_239821/gse239821_data/239821_output"

setwd("~/Documents/nyu/spring/ms_thesis/GSE_239821/gse239821_data")

data <- read_excel(file.path(input_dir, "GSE239821_Raw_counts.xlsx"))

# try to get metadata from geo
gse_meta <- getGEO(GEO = "GSE239821", GSEMatrix = TRUE)

metadata <- pData(phenoData((gse_meta[[1]])))

metadata_updated <- metadata %>% 
  select(1,8, 11) %>% 
  rename(sample = title) %>% 
  rename(treatment = characteristics_ch1.1) %>% 
  rename(cell_source = source_name_ch1) %>% 
  mutate(sample = gsub("Human blastocyst, ", "", sample)) %>% 
  mutate(treatment = gsub("treatment: ", "", treatment)) %>% 
  mutate(sample = gsub("control ([0-9]+)", "CT\\1", sample)) %>%
  mutate(sample = gsub("extracellular vesicles ([0-9]+)", "EV\\1", sample))



count_matrix = as.matrix(data[, -1])
rownames(count_matrix) <- data$`Gene ID`
head(rownames(count_matrix))  # Should display your gene IDs
head(colnames(count_matrix))  # Should display your sample names: 'CT1', 'EV1', etc.

metadata_updated <- metadata_updated %>% 
  arrange(match(sample, colnames(data)[-1]))

all(colnames(count_matrix) == metadata_updated$sample)  # Should return TRUE

rownames(metadata_updated) <- metadata_updated$sample  # This makes sure the row names match the count matrix columns
metadata_updated <- metadata_updated[, -which(names(metadata_updated) == "sample")]  # Remove the 'sample' column

# DESeq2
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("DESeq2")

library(DESeq2)


### gene IDs are specific numbers for genes, which are ensembl number or entrezID,
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
gene_names <- mapIds(org.Hs.eg.db, keys = rownames(count_matrix), 
                     column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
rownames(count_matrix) <- gene_names
count_matrix <- count_matrix[!is.na(rownames(count_matrix)), ]

metadata_updated$treatment <- relevel(factor(metadata_updated$treatment), ref = "None")
## now dds the count_matrix
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = metadata_updated,
                              design = ~ treatment,
                              tidy = FALSE)
dds <- DESeq(dds)
res <- results(dds)

EnhancedVolcano(res,
                lab = rownames(res),
                x = "log2FoldChange",
                y = "padj",
                pCutoff = 0.05,
                FCcutoff = 1,
                title = "EV vs CT (None)",
                subtitle = "Endometrial EV Exposure Effect",
                pointSize = 2.5,
                labSize = 3.0)

# try an MA plot
plotMA(res, ylim = c(-5, 5), main = "MA Plot: EV vs CT")

# MA plot with top 10 genes
plotMA(res, ylim = c(-5, 5), alpha = 0.05, main = "MA Plot: EV vs CT")

# Label top 10 most significant genes
top_genes <- head(order(res$padj), 10)

# Add gene labels
with(res[top_genes, ], {
  text(x = baseMean, y = log2FoldChange, 
       labels = rownames(res)[top_genes],
       pos = 3, cex = 0.7, col = "red")
})

## do I need this?
sig_genes <- outcome[which(outcome$padj < 0.05), ]
head(sig_genes)

# vizz
library(pheatmap)
sig_genes <- rownames(outcome)[which(outcome$padj < 0.05)]
sig_counts <- counts(dds, normalized = TRUE)[sig_genes, ]
pheatmap(log2(sig_counts + 1), cluster_rows = TRUE, cluster_cols = TRUE, 
         show_rownames = TRUE, show_colnames = TRUE, main = "Please work")

# I want to look at specific genes
immune_markers <- c("PTPRC", "CD3E", "CD4", "CD8A", "NCAM1", 
                    "CD28", "FOXP3", "IL2RA", "CD19", 
                    "CD86", "CD14", "ITGAM", "CSF1R", "HOX10", "HOX11", "ITGAV",
                    "IL6", "TNF", "TGFB1", "ERVW-1", "ERVFRD-1", "ERVV-1", "PLAC8", "SERPINE1", "SERPINE2")

marker_counts <- count_matrix[rownames(count_matrix) %in% immune_markers, ]
head(marker_counts)
## 

pheatmap(
  log2(marker_counts + 1),  # Log-transform for better visualization
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,  # Display gene names
  show_colnames = TRUE,  # Display sample names
  main = "Heatmap of Immune Marker Genes"
)

## now that it worked, lets narrow down further ot get the sig genes
# Convert DESeq2 results to data frame
res_df <- as.data.frame(outcome)

# Map your numeric Gene IDs to Gene Symbols
gene_names <- mapIds(
  org.Hs.eg.db,
  keys = rownames(res_df),        # This should be your numeric IDs
  column = "SYMBOL",               # We want the gene symbols
  keytype = "ENTREZID",            # Assuming your IDs are ENTREZ IDs
  multiVals = "first"              # Return the first match if multiple are found
)

res_df$GeneName <- gene_names
sum(!is.na(res_df$GeneName))  # Number of successfully mapped genes
sum(is.na(res_df$GeneName))
# Remove rows with NA GeneNames (if not done already)
res_df_filtered <- res_df[!is.na(res_df$GeneName), ]

immune_markers_2 <- c( 
  "FOXP3", "CD19", 
  "ITGAM", "CSF1R", "HOX10", "HOX11", "ITGAV",
  "IL6", "TGFB1", "ERVW-1", "ERVFRD-1", "ERVV-1", "PLAC8", "SERPINE1", "SERPINE2",
  "CSF1R", "CHN2", "TMEM45B", "GLRX", "DUSP6", "S100A4", "SEC14L4", "SNTB1", "GPR155", "CAV1", "GSTA3",  
  "ACOXL", "MT1H", "DAB2IP", "NUCB2", "CSKMT", "ABHD12B", "LINC00479", "PDK1", "HMOX1")

marker_counts2 <- count_matrix[rownames(count_matrix) %in% immune_markers_2, ]
head(marker_counts2)


pheatmap(
  log2(marker_counts2 + 1),  # Log-transform for better visualization
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,  # Display gene names
  show_colnames = TRUE,  # Display sample names
  main = "Heatmap of Immune Marker Genes"
)

# Convert DESeq2 results to a data frame
res_df <- as.data.frame(outcome)

# Add Gene Names if they are not already present
if (!"GeneName" %in% colnames(res_df)) {
  res_df$GeneName <- rownames(res_df)
}

# Remove rows with NA in the adjusted p-value column (if present)
res_df_filtered <- res_df[!is.na(res_df$padj), ]

# Sort by p-value (or padj if you prefer) to get the most significant genes
top_10_genes <- res_df_filtered %>%
  arrange(padj) %>%
  head(20)

# View the result as a table
print(top_10_genes)

# Install and load the kableExtra package
install.packages("kableExtra")
library(kableExtra)

# Generate a polished table
top_10_genes %>%
  mutate(pvalue = round(as.numeric(pvalue), 4),
         padj = round(as.numeric(padj), 4)) %>%
  kbl(caption = "Top 10 Most Significant Genes from DESeq2 Analysis") %>%
  kable_styling(full_width = F)

top_10_genes <- top_10_genes %>%
  mutate(pvalue = ifelse(as.numeric(pvalue) < 0.05, "<0.05", as.character(pvalue)),
         padj = ifelse(as.numeric(padj) < 0.05, "<0.05", as.character(padj)))

# Generate a polished table with alternate row colors
top_10_genes %>%
  kbl(caption = "Top 10 Most Significant Genes from DESeq2 Analysis") %>%
  kable_styling(full_width = F, bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
  row_spec(0, bold = TRUE, color = "white", background = "#6A0DAD") %>%  # Purple Header
  row_spec(seq(1, 20, 2), background = "#FFF5E1") %>%  # Cream rows (odd)
  row_spec(seq(2, 20, 2), background = "white")        # White rows (even)

# pathway enrichment analysis
library(dplyr)

# Convert res to a dataframe
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Filter by padj and log2FC
up_in_EV <- res_df %>%
  filter(padj < 0.05, log2FoldChange > 1)

up_in_CT <- res_df %>%
  filter(padj < 0.05, log2FoldChange < -1)

library(org.Hs.eg.db)
library(clusterProfiler)

ev_genes <- bitr(up_in_EV$gene, fromType = "SYMBOL",
                 toType = "ENTREZID", OrgDb = org.Hs.eg.db)

ct_genes <- bitr(up_in_CT$gene, fromType = "SYMBOL",
                 toType = "ENTREZID", OrgDb = org.Hs.eg.db)

go_ev <- enrichGO(gene = ev_genes$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP", # Biological Process
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  readable = TRUE)

go_ct <- enrichGO(gene = ct_genes$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  readable = TRUE)
# EV-enriched pathways
barplot(go_ev, showCategory = 10, title = "GO BP: Enriched in EV")

# CT-enriched pathways
barplot(go_ct, showCategory = 10, title = "GO BP: Enriched in CT")

pathways <- go_ev@result %>%
  filter(grepl("implantation|adhesion|matrix|immune|stress|oxidative|cytokine|pluripotency", 
               Description, ignore.case = TRUE)) %>%
  arrange(p.adjust)

dotplot(go_ev, showCategory = pathways$Description)

# LIGAND - RECEPTOR ANALYSIS
"
1. First get the receptors unique to EV and control
2. Perform Differential expression on these receptors
"
if (!requireNamespace("OmnipathR", quietly = TRUE)) {
  BiocManager::install("OmnipathR")
}

# isolate the EV and control from the counts matrix
counts_ev <- count_matrix[, grepl("EV", colnames(count_matrix))]
counts_control <- count_matrix[, grepl("CT", colnames(count_matrix))]

library(OmnipathR)

lr_interactions <- curated_ligand_receptor_interactions()

# Normalize separately
cpm_ev <- cpm(counts_ev)
cpm_control <- cpm(counts_control)

# Get expressed gene names
expressed_genes_ev <- rownames(cpm_ev)[rowSums(cpm_ev > 1) > 0]
expressed_genes_control <- rownames(cpm_control)[rowSums(cpm_control > 1) > 0]

# Define receptor and ligand genes from OmnipathR
receptor_genes <- unique(toupper(lr_interactions$target_genesymbol))
ligand_genes <- unique(toupper(lr_interactions$source_genesymbol))

# Get receptors actually expressed
expressed_receptors_ev <- intersect(expressed_genes_ev, receptor_genes)
expressed_receptors_control <- intersect(expressed_genes_control, receptor_genes)

shared_receptors <- intersect(expressed_receptors_ev, expressed_receptors_control)
ev_only_receptors <- setdiff(expressed_receptors_ev, expressed_receptors_control)
control_only_receptors <- setdiff(expressed_receptors_control, expressed_receptors_ev)

# Deg in shared receptors
dge_matrix <- count_matrix[rownames(count_matrix) %in% shared_receptors, ]

group <- ifelse(grepl("EV", colnames(dge_matrix)), "EV", "Control")
group <- factor(group)

"levels(group)
[1] Control EV

"

# edgeR
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("edgeR")

library(edgeR)

# Create DGEList
dge <- DGEList(counts = dge_matrix, group = group)


# Normalize and estimate
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge)
fit <- exactTest(dge)  # fisher's exact test
 
# Get top DE genes
dge_results <- topTags(fit, n = Inf)$table

library(EnhancedVolcano)

EnhancedVolcano(dge_results,
                lab = rownames(dge_results),
                x = 'logFC',
                y = 'PValue',
                pCutoff = 0.05,
                FCcutoff = 1)


library(pheatmap)

# Select top DE genes
top_genes <- rownames(dge_results)[1:30]  # or based on adj.P.Val < 0.05

# Normalize CPM for plotting
cpm_matrix <- cpm(dge, log = TRUE)
heatmap_data <- cpm_matrix[top_genes, ]

# Annotate groups
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(heatmap_data)

pheatmap(heatmap_data, scale = "row", annotation_col = annotation_col)

# Select top N most variable significant receptors
top_genes <- rownames(head(sig_receptors[order(-abs(sig_receptors$logFC)), ], 20))
mat <- cpm_receptors[top_genes, ]  # or log2(receptor_counts + 1)

# Create annotation for CT vs EV
group <- ifelse(grepl("^CT", colnames(mat)), "CT", "EV")
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(mat)

pheatmap(mat,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main = "Top 20 DE Receptors")

# export ligands as csv 
# Only include interactions where the receptor is in your shared receptor list
filtered_lr <- lr_interactions[lr_interactions$target_genesymbol %in% shared_receptors, ]
ligand_genes <- unique(toupper(filtered_lr$source_genesymbol))

# change to df
ligand_df <- data.frame(Ligand = ligand_genes)

write.csv(ligand_df, "shared_receptor_matched_ligands.csv", row.names = FALSE)

# further analysis
" 
1. group EV and CT separately and check for receptor DGE
2. compare the gene counts deg with the receptor deg"

# de for ev and ct separate
# Assign embryo condition labels
group <- ifelse(grepl("^EV", colnames(receptor_counts)), "EV", "CT")

# edgeR workflow
dge <- DGEList(counts = receptor_counts, group = group)
dge <- calcNormFactors(dge)
design <- model.matrix(~ group)
dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit)
receptor_de <- topTags(lrt, n = Inf)$table

# get teh receptors in EV and control (up in EV and down in control)
sig_receptors_up_in_EV <- subset(sig_receptors, logFC > 0 & FDR < 0.05)
sig_receptors_up_in_control <- subset(sig_receptors, logFC < 0 & FDR < 0.05)

cat(rownames(sig_receptors_up_in_EV), sep = '", "')
cat(rownames(sig_receptors_up_in_control), sep = '", "')
