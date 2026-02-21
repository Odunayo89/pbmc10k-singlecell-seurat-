🧬 Single-Cell RNA-seq Analysis of PBMC (10k Cells) using Seurat
📌 Overview

This project presents a complete single-cell RNA-seq (scRNA-seq) analysis workflow using Seurat in R.
The dataset consists of ~10,000 peripheral blood mononuclear cells (PBMCs) from a healthy donor.

The goal is to identify distinct immune cell populations based on transcriptional profiles.

🔬 Workflow
1. Data Processing

Loaded 10x Genomics HDF5 data

Extracted gene expression matrix (RNA assay only)

2. Quality Control

Filtered cells based on:

Number of detected genes (nFeature_RNA)

Total counts (nCount_RNA)

Mitochondrial gene percentage (%MT)

3. Normalization & Feature Selection

Log-normalization to correct sequencing depth

Identification of 2,000 highly variable genes

4. Dimensionality Reduction

Principal Component Analysis (PCA)

Elbow plot used to select significant PCs

5. Clustering

Graph-based clustering (Louvain algorithm)

Resolution = 0.5

6. Visualization

UMAP for low-dimensional visualization

7. Marker Gene Identification

Differential expression to identify cluster-specific genes

8. Cell Type Annotation

Clusters were annotated based on canonical marker genes.
