library(dplyr)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(readxl)
library(plyr)
library(magrittr)

input_dirs <- list(
  lh = "/gpfs/data/proteomics/home/rs9522/OneDrive_1_14-11-2024/GSE250130/lh",
  pgt = "/gpfs/data/proteomics/home/rs9522/OneDrive_1_14-11-2024/GSE250130/pgt",
  rif = "/gpfs/data/proteomics/home/rs9522/OneDrive_1_14-11-2024/GSE250130/rif"
)

output_dir <- "/gpfs/data/proteomics/home/rs9522/OneDrive_1_14-11-2024/GSE250130/merged_groups"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

setwd("/gpfs/data/proteomics/home/rs9522/OneDrive_1_14-11-2024/GSE250130")

process_directory <- function(input_dir) {
  pat_dirs <- list.dirs(input_dir, recursive = FALSE, full.names = TRUE)
  pat_dirs <- as.list(pat_dirs)
  
  seurat_list <- lapply(pat_dirs, function(dir) {
    print(paste("Processing:", dir))
    
    obj <- (CreateSeuratObject(counts = Read10X(dir), project = basename(dir)))
    
    })
  
}

# reads only the particular directories
seurat_list_lh <- process_directory(input_dirs$lh)
seurat_list_pgt <- process_directory(input_dirs$pgt)
seurat_list_rif <- process_directory(input_dirs$rif)


merged_obj <- function(seurat_list) {
  if (length(seurat_list) > 1) {
    return(merge(seurat_list[[1]], y = seurat_list[-1]))
  } else if (length(seurat_list) == 1) {
    return(seurat_list[[1]])
  } else {
    return(NULL)
  }
}


# calls the function "merged_obj" through seurat_list and created separate obj
lh_obj <- merged_obj(seurat_list_lh)
pgt_obj <- merged_obj(seurat_list_pgt)
rif_obj <- merged_obj(seurat_list_rif)

# merge all three now
obj <- merge(lh_obj, y = c(pgt_obj, rif_obj))

output_path <- file.path(output_dir, "rif_obj_apr22.rds")

# save in the output dir and read in the future
saveRDS(obj, file = output_path)

obj <- readRDS(output_path)

obj[["percent.mt"]] <- PercentageFeatureSet(object = obj, pattern = "^MT-")

# quality metric
obj <- obj[, obj$nFeature_RNA > 250 & obj$nCount_RNA < 100000 & obj$nCount_RNA > 200 & obj$percent.mt < 20]

# normalise 
obj <- obj %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 1000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 30) %>% 
  FindNeighbors(dims = 1:30, k.param = 20) %>% 
  FindClusters(resolution = 1) %>% 
  RunUMAP(dims = 1:30)

# save in the output dir and read in the future
saveRDS(obj, file = output_path)

obj <- readRDS(output_path)


######################################## trying annotation from paper #############################
# intergrate meta data
# I have 4 excel sheets; m4 - patient data; m5 - gene clusters in different phases; m6 - pathway enrichment; m7 deg for RIF patient

library(magrittr)

patient_metadata <- read_excel("/gpfs/data/proteomics/home/rs9522/OneDrive_1_14-11-2024/GSE250130/GSE250130_metadata/41467_2024_55419_MOESM4_ESM.xlsx",
                               sheet = 1) %>% as.data.frame()

head(rownames(obj@meta.data))
head(obj@meta.data)

# both to lowercase
patient_metadata$`Sample ID` <- tolower(patient_metadata$`Sample ID`)
obj$orig.ident <- tolower(obj$orig.ident)

head(colnames(obj@meta.data))
head(rownames(obj@meta.data))

# now to join the metadata with the obj, I am doing all this circus of renaming

# uglyyy namess ughh
pgt_mapping <- c(
  "pat_1-pgt" = "lh7_pgt_1",
  "pat_2-pgt" = "lh7_pgt_2",
  "pat_3-pgt" = "lh7_pgt_3"
)

# rename obj$orig.ident
obj$orig.ident <- mapvalues(
  obj$orig.ident,
  from = names(pgt_mapping),
  to = pgt_mapping
)

unique(obj$orig.ident) # renaming worked

colnames(patient_metadata)[colnames(patient_metadata) == "Sample ID"] <- "orig.ident"

# check colnames 
head(colnames(patient_metadata))

# join metadata
metadata <- obj@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  left_join(patient_metadata, by = "orig.ident") %>%
  tibble::column_to_rownames("barcode")

# Assign back
obj@meta.data <- metadata

# rename donor source col
obj$`Donor source` <- recode(obj$`Donor source`, 
                             "fertile control (public)" = "control",
                             "fertile control (IVF-PGT)" = "IVF", 
                             "IVF patient (RIF)" = "RIF")


# please save - after adding metadata
rif_meta_path <- file.path(output_dir, "obh_june4.rds")
saveRDS(obj, file = rif_meta_path)

obj <- readRDS(rif_meta_path)

# convert colnames
colnames(obj@meta.data) <- make.names(colnames(obj@meta.data))


####################################---Blueprint encode annotations---###############################################

# after reading the obj normalise (this is the metadata intergrated obj)
library(magrittr)
obj <- obj %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 1000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 30) %>% 
  FindNeighbors(dims = 1:30, k.param = 20) %>% 
  FindClusters(resolution = 1) %>%
  RunUMAP(dims = 1:30)

# Rename metadata column
colnames(obj@meta.data)[colnames(obj@meta.data) == "Donor source"] <- "Donor_source"

DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = FALSE)


# cell annotation
library(SingleR)
library(celldex)
library(SingleCellExperiment)

#trying join layers fun
anno_obj <- obj

anno_join <- JoinLayers(anno_obj)

# trying BlueprintENCODE

ref <- celldex::BlueprintEncodeData(ensembl = FALSE)

counts_mat <- LayerData(anno_join, assay = "RNA", layer = "counts")

sce_counts <- SingleCellExperiment(
  assays = list(counts = counts_mat),
  colData = anno_join@meta.data
)

logcounts(sce_counts) <- log1p(counts(sce_counts))

blue_pred <- SingleR(test = sce_counts, ref = ref, labels = ref$label.main)

anno_join$blueprint_label <- blue_pred$labels

# for this particular annotation, I annotated on the non processed obj 
# so now onto our pre processing steps
anno_join[["percent.mt"]] <- PercentageFeatureSet(object = anno_join, pattern = "^MT-")

# quality metric
anno_join <- anno_join[, anno_join$nFeature_RNA > 250 & anno_join$nCount_RNA < 100000 & anno_join$nCount_RNA > 200 & anno_join$percent.mt < 20]

# normalise 
anno_join <- anno_join %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 1000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 30) %>% 
  FindNeighbors(dims = 1:30, k.param = 20) %>% 
  FindClusters(resolution = 1) %>% 
  RunUMAP(dims = 1:30)

DimPlot(anno_join, reduction = "umap", group.by = "blueprint_label", label = TRUE)
FeaturePlot(anno_join, reduction = "umap", features = "EPCAM")
FeaturePlot(anno_join, reduction = "umap", features = "Cd45")

rif_bp <- file.path(output_dir, "anno_bp_jul8.rds")
saveRDS(anno_join, file = rif_bp)

obj<- readRDS(rif_bp)

library(ExperimentHub)
eh <- ExperimentHub()
removeCache(eh)  # clears corrupted cached files

# look for stromal markers
FeaturePlot(obj, features = c("PRL", "FOXO1", "VIM", "COL1A1"))
DimPlot(obj, group.by = "seurat_clusters", label = TRUE)

# Create a new column for modified labels
obj$blueprint_label_modified <- obj$blueprint_label

# Assign 'Stromal cells' to selected clusters
stromal_clusters <- c(1, 2, 5, 11, 12, 25)
obj$blueprint_label_modified[obj$seurat_clusters %in% stromal_clusters] <- "Stromal cells"

# reassigning cell types
obj$blueprint_label_broad <- obj$blueprint_label_modified

obj$blueprint_label_broad[obj$blueprint_label_modified %in% c("Adipocytes", "Chondrocytes", "Mesangial cells", "Skeletal muscle")] <- "Stromal cells"
obj$blueprint_label_broad[obj$blueprint_label_modified %in% c("Astrocytes", "Melanocytes")] <- "Neural"
obj$blueprint_label_broad[obj$blueprint_label_modified == "Erythrocytes"] <- "Blood"
obj$blueprint_label_broad[obj$blueprint_label_modified == "Keratinocytes"] <- "Epithelial"

DimPlot(obj, group.by = "blueprint_label_broad", label = TRUE, repel = TRUE)

# save the reassigned lables 
rif_anno <- file.path(output_dir, "cell_type_jul11.rds")
saveRDS(obj, file = rif_anno)

obj<- readRDS(rif_anno)

# blue print label has a better annotation set than the hpca

# pseudotime 
# Subset for relevant cell types
cell_types_of_interest <- c("Monocytes", "Macrophages", "NK cells")

# Only CONTROL group with full LH+ timepoints
control_cells <- subset(obj,
                        subset = blueprint_label_broad %in% cell_types_of_interest &
                          Donor_source == "control")

colnames(obj@meta.data)[colnames(obj@meta.data) == "Time of biopsy"] <- "Timepoint"

# Only IVF and RIF at LH+7
ivf_rif_cells <- subset(obj,
                        subset = blueprint_label_broad %in% cell_types_of_interest &
                          Donor_source %in% c("IVF", "RIF") &
                          Timepoint == "LH+7")

library(monocle3)
library(SeuratWrappers)

# conver subset obj to cds
cds_control <- as.cell_data_set(control_cells)
cds_control <- cluster_cells(cds_control)
cds_control <- learn_graph(cds_control)

colnames(colData(cds_control))[colnames(colData(cds_control)) == "Time.of.biopsy"] <- "Timepoint"

# Visualize to pick root if needed:
plot_cells(cds_control, color_cells_by = "Timepoint")

# Use LH+3 as root
# Extract LH+3 cell barcodes
lh3_cells <- rownames(subset(colData(cds_control), Timepoint == "LH+3"))
cds_control <- order_cells(cds_control, root_cells = lh3_cells)

plot_cells(cds_control, color_cells_by = "pseudotime", label_groups_by_cluster = FALSE)


cds_control@colData$pseudotime <- pseudotime(cds_control)
boxplot(pseudotime ~ Timepoint, data = as.data.frame(colData(cds_control)))

plot_cells(cds_control, color_cells_by = "blueprint_label_broad") #original time point

# psudotime for the ivf and rif group
cds_ivf_rif <- as.cell_data_set(ivf_rif_cells)
cds_ivf_rif <- cluster_cells(cds_ivf_rif)
cds_ivf_rif <- learn_graph(cds_ivf_rif)

# Pick random root for now since no prior timepoint
some_cells <- colnames(cds_ivf_rif)[1:20]
cds_ivf_rif <- order_cells(cds_ivf_rif, root_cells = some_cells)

# Visualize
plot_cells(cds_ivf_rif, color_cells_by = "pseudotime", label_groups_by_cluster = FALSE)

plot_cells(cds_ivf_rif, color_cells_by = "blueprint_label_broad")

########################CCC using LIANA ######################

devtools::install_github("saezlab/liana")
# Install latest LIANA with Seurat support
devtools::install_github("saezlab/liana", ref = "dev")

library(OmnipathR)

# get the ligands for the rif dataset; we got the receptor expression for embryo
ligrec <- curated_ligand_receptor_interactions()
ligand_genes <- unique(toupper(ligrec$source_genesymbol))


# Ensure gene names are uppercase for matching
rownames(obj) <- toupper(rownames(obj))

# Subset for ligands
ligand_genes_in_obj <- intersect(ligand_genes, rownames(obj))
length(ligand_genes_in_obj)  # how many ligands expressed in your data

# Optionally, create a module score
obj <- AddModuleScore(obj, features = list(ligand_genes_in_obj), name = "LigandScore")

DotPlot(obj, features = ligand_genes_in_obj[1:20], group.by = "blueprint_label_broad", cols = c("lightgrey", "darkred")) +
  RotatedAxis() +
  labs(title = "Top Ligand Expression by Cell Type")

# by condition
Idents(obj) <- obj$Donor_source
ligand_de <- FindMarkers(obj, ident.1 = "RIF", ident.2 = "control",
                         features = ligand_genes_in_obj,
                         logfc.threshold = 0.25, min.pct = 0.1)

# View top results
head(ligand_de)

# Add gene names as a column
ligand_de$gene <- rownames(ligand_de)

# Filter significant DE ligands (FDR < 0.05)
sig_ligands <- subset(ligand_de, p_val_adj < 0.05)

ligand_de$significance <- ifelse(ligand_de$p_val_adj < 0.05, "FDR < 0.05", "Not Sig")

top_genes <- head(sig_ligands$gene, 5)

DotPlot(obj, features = top_genes, group.by = "Donor_source") +
  RotatedAxis() +
  labs(title = "Top DE Ligands by Donor Source")

VlnPlot(obj, features = top_genes, group.by = "Donor_source", pt.size = 0.1)

top_ligands <- head(sig_ligands[order(-abs(sig_ligands$avg_log2FC)), ], 20)

ggplot(top_ligands, aes(x = reorder(gene, avg_log2FC), y = avg_log2FC, fill = avg_log2FC > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("firebrick", "steelblue"), labels = c("Down in RIF", "Up in RIF")) +
  labs(title = "Differentially expressed ligands", x = "Gene", y = "log2 Fold Change") +
  theme_minimal()

# getting the receptor de form the embryos
# Get your DE ligands (from RIF endometrium)
ligands <- rownames(sig_ligands)  # from your barplot

# Get DE receptors (from embryo dataset heatmap or DE result)
receptors <- c("TNFSF8", "APCDD1", "ADGRB2", "GREM2", "ALK", "CAV1", "CSF1R", "FLT1", 
               "IL21R", "ITGA8", "ACVR2B", "NCR1", "CD82", "TRAF2", "TRADD", 
               "ADGRB1", "AQP6", "VIPR2", "RXFP4", "TSPAN17")

# Filter OmniPath pairs that match
lr_pairs_filtered <- subset(ligrec,
                            toupper(source_genesymbol) %in% toupper(ligands) &
                              toupper(target_genesymbol) %in% toupper(receptors))

lr_pairs_filtered[, c("source_genesymbol", "target_genesymbol")]

# vizualise the interaction
library(igraph)

g <- graph_from_data_frame(lr_pairs_filtered[, c("source_genesymbol", "target_genesymbol")], directed = TRUE)

V(g)$type <- ifelse(V(g)$name %in% ligands, "Ligand", "Receptor")
V(g)$color <- ifelse(V(g)$type == "Ligand", "tomato", "skyblue")

plot(g,
     vertex.label.cex = 0.8,
     edge.arrow.size = 0.5,
     layout = layout_with_fr(g),
     main = "Ligand–Receptor Interactions: RIF ↔ Embryo")

# make circos plot 
install.packages("circlize")
library(circlize)


# Example input: chord_data with Ligand–Receptor columns
circos.clear()
chordDiagram(chord_data,
             grid.col = grid.col,               # if you're using color mapping
             transparency = 0.5,
             directional = 1,
             direction.type = "arrows",
             annotationTrack = "grid",          # ensures labels show
             preAllocateTracks = list(track.height = 0.1))  # ensures space for labels

# Add labels manually if still missing:
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], 
              CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.6)
}, bg.border = NA)

#### L- R for each of the groups in the dataset
"Ligand source (Endometrium)  Receptor source (Embryo)

RIF - Receptors up in EV vs control


IVF - Receptors up in EV vs control

control - Receptors up in EV vs control

RIF - Receptors up in control vs EV

IVF  - Receptors up in control vs EV

control - Receptors up in control vs EV
"

ev_receptors <- c("CSF1R", "TSPAN17", "CAV1", "TRADD", "AQP6", "TRAF2", 
                 "VLDLR", "FLT1", "NCR1", "GREM2", "RORA", "SCN8A", "ACVR2B", 
                 "CD82", "SIGIRR", "TRPV6", "CAP1", "ITGB1", "VIPR2", "ACTR2", 
                 "CXCR6", "ADGRB2", "RXFP4", "ALK", "AVPR1A", "ITGA8", "EGFR", 
                 "IL21R", "CKAP4", "ADGRB1", "CXCR4", "FZD4", "LPAR3", "ENG", 
                 "RXFP1", "STRA6", "ADIPOR1", "HBEGF", "AXL", "TNFRSF12A", "NECTIN3")


control_receptors <- c("APCDD1", "KLRK1", "SLC4A11", "S1PR3", "CD36", "EDA2R", "F2R", 
"LIFR", "TNFSF8", "ADGRE5", "TNFRSF10C", "SIGLEC5", "HAVCR2", "CFTR", "ADAM11", "MYLK2")

rif_obj <- subset(obj, subset = Donor_source == "RIF")

# Match RIF ligands to EV-upregulated embryo receptors
lr_pairs <- curated_ligand_receptor_interactions()

# Filter ligand–receptor pairs
matched_ev <- subset(lr_pairs,
                     toupper(source_genesymbol) %in% ligand_genes_in_obj &
                       toupper(target_genesymbol) %in% ev_receptors
)

matched_ctrl <- subset(lr_pairs,
                       toupper(source_genesymbol) %in% ligand_genes_in_obj &
                         toupper(target_genesymbol) %in% control_receptors
)

########################
library(SingleCellExperiment)
library(liana)

obj$celltype <- obj$blueprint_label_broad

# Stromal and immune
stromal_immune_obj <- subset(obj, subset = celltype %in% c("Stromal cells", "Macrophages", "Monocytes", "NK cells"))

# convert obj to sce
sce_obj <- as.SingleCellExperiment(stromal_immune_obj)

colData(sce_obj)$celltype <- stromal_immune_obj$blueprint_label_broad
colData(sce_obj)$Donor_source <- stromal_immune_obj$Donor_source

liana_result <- liana_wrap(
  sce = sce_obj,
  group_by = "Donor_source",     # 'control', 'IVF', 'RIF'
  assay = "logcounts",           # use log-normalized expression
  idents_col = "celltype",       # cell type column you just added
  method = "consensus"           # aggregate across CCC methods
)
head(liana_result)

stromal_to_immune <- liana_result %>%
  filter(source == "Stromal cells" & target %in% c("Macrophages", "Monocytes", "NK cells"))

library(ggplot2)

ggplot(stromal_to_immune, aes(x = ligand, y = receptor, size = magnitude, color = target)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

liana_ctrl <- liana_wrap(
  seurat_obj = subset(stromal_immune_obj, subset = Donor_source == "control"),
  group_by = NULL,
  assay = "RNA",
  slot = "data",
  celltype_col = "celltype"
)

liana_rif <- liana_wrap(
  seurat_obj = subset(stromal_immune_obj, subset = Donor_source == "RIF"),
  group_by = NULL,
  assay = "RNA",
  slot = "data",
  celltype_col = "celltype"
)

# Compare ligand–receptor scores
inner_join(liana_ctrl, liana_rif, by = c("ligand", "receptor", "source", "target"), suffix = c("_ctrl", "_rif")) %>%
  mutate(diff = magnitude_rif - magnitude_ctrl)

###################################################################

# immune subset from obj
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

obj$cell_type <- "Unknown"

obj <- DietSeurat(obj, assays = "RNA", layers = "data") 

detected_genes <- rownames(GetAssayData(obj, assay = "RNA", layer = "data"))

# now, assign cell types to their markers
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

immune_cells <- subset(obj, subset = cell_type != "Unknown")
unique(immune_cells$cell_type)

#to increase size for parallel processing
options(future.globals.maxSize = 4 * 1024^3)  # Set to 4GB

#renormalize
immune_obj <- immune_cells %>%
  SCTransform(vars.to.regress = c("percent.mt"), verbose = FALSE) %>%
  RunPCA(npcs = 25) %>%
  RunUMAP(dims = 1:25) %>%
  FindNeighbors(dims = 1:25, k.param = 20) %>%
  FindClusters(resolution = 0.5)

# Select top 3000 variable genes
VariableFeatures(immune_obj) <- head(VariableFeatures(immune_obj), 3000)

#vizualise these cuties
DimPlot(immune_obj, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE)

# FeaturePlot(immune_obj, features = c("PTPRC", "CD4", "CD8A", "ITGAX", "MS4A2", "MS4A1"))

p1 <- FeaturePlot(immune_obj, features = "PTPRC", cols = c("gray", "red"))
p2 <- FeaturePlot(immune_obj, features = "CD4", cols = c("gray", "blue"))
p3 <- FeaturePlot(immune_obj, features = "CD8A", cols = c("gray", "green"))
p4 <- FeaturePlot(immune_obj, features = "ITGAX", cols = c("gray", "purple"))
p5 <- FeaturePlot(immune_obj, features = "MS4A2", cols = c("gray", "orange"))
p6 <- FeaturePlot(immune_obj, features = "MS4A1", cols = c("gray", "brown"))

# Arrange in a 2-row grid - chatgpt for the win
(p1 | p2 | p3) / (p4 | p5 | p6)




