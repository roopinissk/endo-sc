# files 
library(dplyr)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(patchwork)


input_dir <- "/gpfs/data/proteomics/home/rs9522/OneDrive_1_14-11-2024/GSE214411"
output_dir <- "/gpfs/data/proteomics/home/rs9522/OneDrive_1_14-11-2024/GSE214411/endo"

dir.create(output_dir, recursive = TRUE)

setwd("/gpfs/data/proteomics/home/rs9522/OneDrive_1_14-11-2024/GSE214411")

rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)

#naming the files
for (file in rds_files) {
  new_name <- gsub("^GSE214411_", "", basename(file))
  file.rename(file, file.path(input_dir, new_name))
}

list.files(path = input_dir, pattern = "\\.rds$", full.names = TRUE) #this works, the file names are eg EMS1_obj

# list the files again since you renamed
rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)

# read rds files
obj_list <- lapply(rds_files, readRDS)

# merge all so
obj <- Reduce(function(x,y) merge(x, y), obj_list)

# Filter cells based on quality criteria
obj <- obj[, obj$nFeature_RNA > 250 & obj$nCount_RNA < 100000 & obj$percent.mt < 20 & obj$nCount_RNA > 200]

# Preprocessing pipeline
obj <- obj %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 25) %>%
  RunUMAP(dims = 1:25) %>%
  FindNeighbors(k.param = 20) %>%
  FindClusters(resolution = 1)

# Save the processed object
saveRDS(obj, file = "merged_seurat_obj.rds")

#if you wanna do other differential exp or something else, simple read the processed file
obj <- readRDS("merged_seurat_obj.rds")

#vizualise
DimPlot(obj, reduction = "pca", group.by = "seurat_clusters", label = TRUE, pt.size = 0.5) +ggtitle("PCA clusters")
ggsave("pca_gse214411", pca_plot)

DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, pt.size = 0.5) +ggtitle("UMAP clusters")
ggsave(umap_plot, "umap_gse214411")

DimPlot(obj, reduction = "umap", group.by = "Lab_ID", label = TRUE, pt.size = 0.5) +ggtitle("UMAP clusters")
ggsave(umap_plot, "umap_gse214411")


#annotating cells for feature plots
#immune cells from paper
markers <- list( 
  Th1 = c("CD3E", "CD4", "IFNG"), #cd3d otherwise
  Th2 = c("CD3E","CD4", "GATA3"),
  Th17 = c("CD3E","CD4", "RORC"),
  Treg = c("CD3E","CD4", "FOXP3"),
  Tc = c("CD3E","CD8A"),
  NK = c("CD3E","NCAM1"),
  Macrophage_M1 = c("CD68", "CD86"),
  Macrophage_M2 = c("CD68", "MRC1"),
  Dendritic = c("ITGAX"),
  Mast = c("MS4A2"),
  B_cells = c("MS4A1") 
  ) 

'1. Add module score
2. expand more genes for anything wiht one'
#add it to metadata; bleh
obj$cell_type <- "Unknown"

#now, assign cell types to their markers
for (cell_type in names(markers)) {
  gene_set <- markers[[cell_type]]
  detected_genes <- rownames(obj[["RNA"]]@data)
  valid_markers <- gene_set[gene_set %in% detected_genes]
  
  if (length(valid_markers) > 0) {
    expr_data <- FetchData(obj, vars = valid_markers)
    expr <- rowSums(expr_data >0.2) >0
    obj$cell_type[expr] <- cell_type
  }
}

obj$cell_type <- factor(obj$cell_type)

#################################################### Blueprint label - Annotations ########################
library(celldex)
library(SingleR)
library(SingleCellExperiment)
ref <- celldex::BlueprintEncodeData()

sce <- as.SingleCellExperiment(obj)


pred <- SingleR(test = sce, ref = ref, labels = ref$label.fine,
                assay.type.test = "logcounts")

obj$Blueprint_main <- pred$pruned.labels

unique(obj$Blueprint_main)


DimPlot(obj, reduction = "umap", group.by = "Blueprint_main", label = TRUE, repel = TRUE)

saveRDS(obj, file = "anno_bp_214411_jul11.rds")
obj <- readRDS("anno_bp_214411_jul11.rds")


########################################


# monaco immune - checking to get more refined annotations using monaco immune
sce_obj <- immune_obj
sce_immune<-as.SingleCellExperiment(sce_obj)

ref_cells_m <- celldex::MonacoImmuneData()
pred_m <- SingleR(test = sce_immune, ref = ref_cells_m, labels = colData(ref_cells_m)$label.main)

# store predicted labels
immune_obj$SingleR_labels_monaco <- factor(pred_m$labels)
immune_obj<-SetIdent(immune_obj, value = "SingleR_labels")
DimPlot(immune_obj,group.by = "SingleR_labels_monaco", reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("MonacoImmuneData")


#change condition to c("mild","moderate","healthy") ## 1-5, 6-15, controls
immune_obj$ASRM_stage<-"mild"
immune_obj$ASRM_stage[immune_obj$Stage=="NA"]<-"healthy"
immune_obj$ASRM_stage[immune_obj$Stage=="6"]<-"moderate"
immune_obj$ASRM_stage<-factor(immune_obj$ASRM_stage)

immune_obj$ASRM_stage_combo<-"mild_moderate"
immune_obj$ASRM_stage_combo[immune_obj$Stage=="NA"]<-"healthy"
immune_obj$ASRM_stage_combo<-factor(immune_obj$ASRM_stage_combo)

Idents(immune_obj)<-"Tumor_histology"
hey_there<-FindMarkers(immune_obj, subset.ident = "Early_proliferative_phase", 
            group.by = "ASRM_stage_combo", assay="RNA",
            ident.1 = "mild_moderate")

View(hey_there)
DimPlot(immune_obj, cells.highlight = WhichCells(immune_obj, expression = PTPRC >=1.25))

options(future.globals.maxSize = 2 * 1024^3) # 2 GB

ptrpc <- subset(immune_obj, cells=WhichCells(immune_obj, expression = PTPRC >=1.25))
DimPlot(ptrpc, reduction='umap', group.by="seurat_clusters")
ptrpc<- ptrpc %>% FindNeighbors() %>% FindClusters()

ptrpc <- ptrpc %>%
  SCTransform(vars.to.regress = c("percent.mt"), verbose = FALSE) %>%
  RunPCA(npcs = 20) %>%
  RunUMAP(dims = 1:20) %>%
  FindNeighbors(dims = 1:20, k.param = 20) %>%
  FindClusters(resolution = 0.5)

DimPlot(ptrpc, group.by = "seurat_clusters")
p1 <- FeaturePlot(ptrpc, features = "CD3E", cols = c("gray", "red"))
p2 <- FeaturePlot(ptrpc, features = "CD4", cols = c("gray", "blue"), order=T)
p3 <- FeaturePlot(ptrpc, features = "CD8A", cols = c("gray", "green"))
p4 <- FeaturePlot(ptrpc, features = "NCAM1", cols = c("gray", "purple"))
p5 <- FeaturePlot(ptrpc, features = "MS4A2", cols = c("gray", "orange"))
p6 <- FeaturePlot(ptrpc, features = "MS4A1", cols = c("gray", "brown"))

# Arrange in a 2-row grid - chatgpt for the win
(p1 | p2 | p3) / (p4 | p5 | p6)

###########################################################

#cd3e+ -- annotate; 
# maybe cd4, cd8, ncam1
# subset just these cells of interest

new_cells <- subset(immune_obj, subset = CD3E > 0 | CD4 > 0 | CD8A > 0 | NCAM1 > 0)
new_cells <- new_cells %>% 
  SCTransform(vars.to.regress = "percent.mt", verbose = FALSE) %>% 
  RunPCA(npcs = 10) %>% 
  RunUMAP(dims = 1:10) %>% 
  FindNeighbors(dims = 1:10, k.param = 20) %>% 
  FindClusters(resolution = 0.5)

DimPlot(new_cells, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
DimPlot(new_cells, reduction = "pca", group.by = "seurat_clusters", label = TRUE, repel = TRUE)

VlnPlot(new_cells, features = c("CD3E", "CD4", "CD8A", "NCAM1", "ITGAX"), group.by = "seurat_clusters", pt.size = 0)

sub_new <- subset(new_cells, cells = WhichCells(new_cells, expression = CD3E > 0.5 |CD4 > 0.5
                                             | CD8A > 0.5 | NCAM1 >0.5 | ITGAX > 0.5))


# Vln plots
p1 <- VlnPlot(sub_new, features = "CD3E", group.by = "seurat_clusters", pt.size = 0)
p2 <- VlnPlot(sub_new, features = "CD4", group.by = "seurat_clusters", pt.size = 0)
p3 <- VlnPlot(sub_new, features = "CD8A", group.by = "seurat_clusters", pt.size = 0)
p4 <- VlnPlot(sub_new, features = "NCAM1", group.by = "seurat_clusters", pt.size = 0)
p5 <- VlnPlot(sub_new, features = "ITGAX", group.by = "seurat_clusters", pt.size = 0)

# Arrange vertically
p1 / p2 / p3 / p4 / p5 + plot_layout(guides = "collect") 

#########################################################

############################ tumor histology ###########################
# early proliferative; mid secretory, late secretory

immune_ems <- subset(immune_obj, subset = ASRM_stage_combo == "mild_moderate")
immune_n <- subset(immune_obj, subset = ASRM_stage_combo == "healthy")

# pre process - EMS
immune_ems <- immune_ems %>% 
  SCTransform(vars.to.regress = "percent.mt", verbose = FALSE) %>% 
  RunPCA(npcs = 10) %>% 
  RunUMAP(dims = 1:10) %>% 
  FindNeighbors(dims = 1:10, k.param = 20) %>% 
  FindClusters(resolution = 0.5)

# pre process - N
immune_n <- immune_n %>% 
  SCTransform(vars.to.regress = "percent.mt", verbose = FALSE) %>% 
  RunPCA(npcs = 10) %>% 
  RunUMAP(dims = 1:10) %>% 
  FindNeighbors(dims = 1:10, k.param = 20) %>% 
  FindClusters(resolution = 0.5)

# let see them
DimPlot(immune_ems, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("EMS seurat")
DimPlot(immune_n, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("N seurat")

# across the phases
ems_early <- subset(immune_ems, subset = Tumor_histology == "Early_proliferative_phase")
ems_mid <- subset(immune_ems, subset = Tumor_histology == "Mid_secretory_phase")
ems_late <- subset(immune_ems, subset = Tumor_histology == "Late_secretory_phase") #not present, idk why

# n
n_early <- subset(immune_n, subset = Tumor_histology == "Early_proliferative_phase")
n_mid <-  subset(immune_n, subset = Tumor_histology == "Mid_secretory_phase")
n_late <- subset(immune_n, subset = Tumor_histology == "Late_secretory_phase")

# plotsss
p1 <- VlnPlot(ems_early, features = c("CD3E", "CD4", "CD8A", "NCAM1", "ITGAX"), 
              group.by = "Tumor_histology", pt.size = 0) + ggtitle("EMS - Early Proliferative Phase")
p2 <- VlnPlot(ems_mid, features = c("CD3E", "CD4", "CD8A", "NCAM1", "ITGAX"), 
              group.by = "Tumor_histology", pt.size = 0) + ggtitle("EMS - Mid_secretory_phase")
p3 <- VlnPlot(ems_late, features = c("CD3E", "CD4", "CD8A", "NCAM1", "ITGAX"), 
              group.by = "Tumor_histology", pt.size = 0) + ggtitle("EMS - Late_secretory_phase")


p1 | p2 + plot_layout(guides = "collect")







#######################################################
obj2<-obj
obj2<-AddModuleScore(obj2,features =  markers, name = names(markers))
rm(obj2)
obj$Stage.points<-obj$Stage
obj$Stage[obj$Stage.points=="NA"]<-"Control"
obj$Stage[obj$Stage.points=="1" | obj$Stage.points=="2"|obj$Stage.points=="3"|obj$Stage.points=="4"|obj$Stage.points=="5"]<-"Mild"
obj$Stage[obj$Stage.points=="6"]<-"Moderate"

obj <- CellCycleScoring(obj, s.features = cc.genes.updated.2019$s.genes, 
                         g2m.features = cc.genes.updated.2019$g2m.genes, assay = "RNA")
DimPlot(obj, group.by = "Stage")

#which cells
#reduce clusters
#immune cells based on staging
#barplot or heatmap -- cell type


#filter unknown cells
obj_filter <- subset(obj, subset = cell_type != "Unknown")

VlnPlot(obj, features = c("CD4", "CD8A", "ITGAX", "MS4A2", "MS4A1"), group.by = "cell_type")

Idents(obj_filter) <- "cell_type"
FindMarkers(obj, ident.1 = "Unknown", only.pos = TRUE)  #seems like mostly fibroblast or stromal cells - chatgpt told me

DimPlot(obj, reduction = "umap", group.by = "cell_type", label = TRUE) +
  ggtitle("Annotated immune cells - UMAP")

Idents(obj) <- "cell_type"  # Set cell type as the identity
immune_DE <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(immune_DE)


#tumor histology DE
Idents(obj) <- "Tumor_histology"
early_vs_late <- FindMarkers(obj, ident.1 = "Early_proliferative_phase", ident.2 = "Late_secretory_phase",
                             min.pct = 0.25, logfc.threshold = 0.25)
head(early_vs_late)


early_vs_late$gene <- rownames(early_vs_late)  


FeaturePlot(obj, features = c("IGFBP5", "SFRP4", "COL1A1", "SH3BGRL3", "BTG1"), reduction = "umap")

#trying to set genes to tumor histology
Idents(obj) <- "Tumor_histology"  # Set menstrual phases as identity

immune_DE <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(immune_DE)

DefaultAssay(obj) <- "RNA"

all_genes <- unique(c(top_markers$gene, rownames(obj)))  # Ensure missing genes are included
obj <- ScaleData(obj, features = all_genes)

heatmap <- DoHeatmap(obj, features = top_markers$gene, assay = "RNA", slot = "counts")
print(heatmap)

immune_genes <- c("IL6", "IL1B", "TNF", "CCL2", "CCL4", "CXCL8", "CXCL12",
                  "CD68", "CD3E", "FOXP3", "NCAM1", "MMP9", "TIMP1",
                  "CD4", "CD8A", "ITGAX")

implantation_genes <- c("WNT4", "ITGAV", "ITGB3", "LIF",
                        "HOXA10", "IL10", "TGFB1")

# Combine immune and implantation genes
selected_genes <- c(immune_genes, implantation_genes)

# Get average expression for these genes across cell types
avg_expression <- AverageExpression(immune_obj, features = selected_genes, group.by = "cell_type")$RNA

# Convert to long format for plotting
dotplot_data <- avg_expression %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "gene") %>%
  pivot_longer(cols = -gene, names_to = "cell_type", values_to = "expression")

menstrual_phase <- as_tibble(menstrual_phase)
colnames(menstrual_phase)

# get DEG from the mid secretory phase
DotPlot(mid_secretory_obj, 
        features = all_genes, 
        group.by = "ASRM_stage_combo") +
  ggtitle("Immune and Implantation Gene Expression in Healthy vs. Mild-Moderate") +
  scale_color_gradient(low = "#FF4500", high = "darkred") +  # Gradient from Orange to Dark Red
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),  # Reduce text size on y-axis for compression
    plot.margin = margin(5, 5, 5, 5),
    panel.spacing.y = unit(0.05, "lines")  # Reduce space between mild_moderate and healthy
  )

# Merge with DE p-values - Dont use
dotplot_data <- dotplot_data %>%
  left_join(dplyr::select(menstrual_phase, gene, p_val), by = "gene") %>%
  mutate(significance = -log10(p_val))  # Convert p-value to -log10 scale

ggplot(dotplot_data, aes(x = cell_type, y = gene)) +
  geom_point(aes(size = significance, color = expression)) +
  scale_color_gradientn(colors = c("#FF6600", "#FF9966", "#FF99CC", "#FF66FF")) + 
  scale_size_continuous(range = c(1, 8)) +  # Adjust dot sizes
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Expression of Immune and Implantation Genes in Immune Cells",
       x = "Immune Cell Type", y = "Gene",
       color = "Log2 Expression", size = "-log10(p-value)")

# Pathway enrichment analysis
library(dplyr)
library(tibble)

# Convert row names to a column
menstrual_phase <- menstrual_phase %>%
  rownames_to_column(var = "gene")

# Select significant genes
sig_genes <- menstrual_phase %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5) %>%
  arrange(desc(avg_log2FC))

# Extract upregulated and downregulated genes
upregulated_genes <- sig_genes %>% filter(avg_log2FC > 0) %>% pull(gene)
downregulated_genes <- sig_genes %>% filter(avg_log2FC < 0) %>% pull(gene)

# Convert gene symbols to Entrez IDs
gene_entrez_ids <- bitr(sig_genes$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # very important

# GO Enrichment for Biological Processes (BP)
go_enrich <- enrichGO(
  gene = gene_entrez_ids$ENTREZID, 
  OrgDb = org.Hs.eg.db, 
  keyType = "ENTREZID", 
  ont = "BP",  # Biological Process----> can also use "MF" for Molecular Function, "CC" for Cellular Component; this is most relevant for tjis dataset
  pAdjustMethod = "BH", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05
)

# Visualize GO results
dotplot(go_enrich, showCategory = 10, title = "GO Enrichment Analysis - Biological Processes")

 