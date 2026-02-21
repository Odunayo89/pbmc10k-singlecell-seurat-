# Single-Cell RNA-seq Analysis of PBMC (10k Cells) using Seurat
## Overview

This project demonstrates a complete single-cell RNA-seq (scRNA-seq) analysis workflow using Seurat in R. The dataset consists of ~10,000 peripheral blood mononuclear cells (PBMCs) from a healthy donor.

The objective is to identify and characterize distinct immune cell populations based on their transcriptional profiles.

### Workflow
1. Data Processing

Loaded 10x Genomics HDF5 dataset

Extracted gene expression matrix (RNA assay)

2. Quality Control

Filtered low-quality cells based on:

Number of detected genes (nFeature_RNA)

Total RNA counts (nCount_RNA)

Mitochondrial gene percentage (%MT)

3. Normalization & Feature Selection

Log-normalization to correct for sequencing depth differences

Identification of 2,000 highly variable genes

4. Dimensionality Reduction

Principal Component Analysis (PCA)

Elbow plot used to determine informative components

5. Clustering

Graph-based clustering (Louvain algorithm)

Resolution parameter set to 0.5

6. Visualization

UMAP used to visualize transcriptional similarity between cells

7. Marker Gene Identification

Differential expression analysis performed to identify cluster-specific genes

8. Cell Type Annotation

Clusters annotated using canonical immune cell marker genes

### Results
UMAP Visualization

Marker Gene Expression (DotPlot)

Cluster Marker Heatmap




### Identified Cell Types
Cluster	Cell Type
0	CD14+ Monocytes
1	CD4 T cells
2	CD8 T cells
3	Inflammatory Monocytes
4	Activated T cells
5	B cells
6	Mixed/Unknown
7	NK cells
8	Rare population
9	CD16+ Monocytes
10	Rare population

### Tools & Libraries

R

Seurat

dplyr

### Dataset

Public dataset from 10x Genomics:
PBMC from a Healthy Donor (~10k cells)

Note: Raw data files are not included. Please download from 10x Genomics and place the .h5 file in the project directory.

### Key Takeaways

Successfully identified major immune cell populations from scRNA-seq data

Demonstrated a full end-to-end Seurat workflow

Translated transcriptional patterns into biological cell identities

👤 Author

Adekunle Ajiboye
