🧬 Single-cell RNA-seq Analysis of PBMC 10k (Seurat Workflow)
Overview

This project demonstrates a complete single-cell RNA-seq analysis workflow using Seurat in R, applied to a publicly available 10x Genomics PBMC dataset.

The analysis includes quality control, normalization, dimensionality reduction, clustering, and biological annotation of immune cell populations.

Workflow

Data loading (10x Genomics HDF5 format)

Quality control (mitochondrial filtering, feature counts)

Normalization and identification of highly variable genes

Dimensionality reduction (PCA)

Graph-based clustering

Visualization using UMAP

Marker gene identification and cluster annotation

Key Results

Identification of major immune cell populations:

CD4 T cells

CD8 T cells

B cells

NK cells

Monocytes

Dendritic cells

Tools Used

R

Seurat

dplyr

Dataset

Publicly available dataset from 10x Genomics:
PBMC from a healthy donor (10k cells)

Author

Adekunle Ajiboye
