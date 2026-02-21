library(Seurat)
library(dplyr)

# 0) Load data (Multiome file -> take Gene Expression only)
counts <- Read10X_h5("filtered_feature_bc_matrix.h5")
names(counts)  # should include "Gene Expression" (and maybe "Peaks")

pbmc10k <- CreateSeuratObject(
  counts = counts[["Gene Expression"]],
  project = "pbmc10k",
  min.cells = 3,
  min.features = 200
)

# 1) QC (recommended)
pbmc10k[["percent.mt"]] <- PercentageFeatureSet(pbmc10k, pattern = "^MT-")
VlnPlot(pbmc10k, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(pbmc10k, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(pbmc10k, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Optional: filter (adjust after inspecting violins/scatters)
pbmc10k <- subset(pbmc10k,
                  subset = nFeature_RNA > 300 &
                    nFeature_RNA < 6000 &
                    percent.mt < 10)

# 2) Normalize + Variable Features
pbmc10k <- NormalizeData(pbmc10k)
pbmc10k <- FindVariableFeatures(pbmc10k, selection.method = "vst", nfeatures = 2000)

# sanity check
length(VariableFeatures(pbmc10k))
head(VariableFeatures(pbmc10k))

# 3) Scale + PCA
pbmc10k <- ScaleData(pbmc10k)
pbmc10k <- RunPCA(pbmc10k)

ElbowPlot(pbmc10k)

# 4) Neighbors + Clusters + UMAP
pbmc10k <- FindNeighbors(pbmc10k, dims = 1:20)
pbmc10k <- FindClusters(pbmc10k, resolution = 0.5)

pbmc10k <- RunUMAP(pbmc10k, dims = 1:20)
DimPlot(pbmc10k, label = TRUE)

# 5) Markers
markers10k <- FindAllMarkers(pbmc10k, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(markers10k)

# Top 10 markers per cluster (for heatmap)
top10 <- markers10k %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

# Ensure these genes are scaled (DoHeatmap uses scale.data)
pbmc10k <- ScaleData(pbmc10k, features = unique(top10$gene))
DoHeatmap(pbmc10k, features = unique(top10$gene)) + NoLegend()

# 6) Quick cell-type marker visualization
FeaturePlot(pbmc10k, features = c(
  "CD3D", "IL7R",        # T cells (IL7R = CD4 T)
  "NKG7", "GNLY",        # NK / cytotoxic
  "MS4A1", "CD79A",      # B cells
  "LYZ", "S100A8",       # Monocytes (S100A8 = classical)
  "FCGR3A", "MS4A7",     # CD16 mono
  "FCER1A", "CST3",      # Dendritic
  "PPBP"                 # Platelets
), ncol = 4)

# DotPlot (clean summary for GitHub)
DotPlot(pbmc10k, features = c(
  "CD3D","IL7R","NKG7","GNLY","MS4A1","CD79A","LYZ","FCGR3A","FCER1A","PPBP"
)) + RotatedAxis()

# 7) (Optional) Rename clusters AFTER you inspect markers
# Fill in based on what lights up in FeaturePlot/DotPlot:
# new.ids <- c("0"="CD14+ Mono", "1"="CD4 T", "2"="NK", ...)
# pbmc10k <- RenameIdents(pbmc10k, new.ids)
# DimPlot(pbmc10k, label = TRUE) + NoLegend()

# 8) Save outputs (recommended for GitHub/reproducibility)
saveRDS(pbmc10k, "pbmc10k_seurat.rds")
write.csv(markers10k, "pbmc10k_markers.csv", row.names = FALSE)