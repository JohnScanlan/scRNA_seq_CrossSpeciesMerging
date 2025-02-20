# scRNA_seq_CrossSpeciesMerging
Combining Human and Mouse scRNA-seq data

This project demonstrates how to integrate single-cell RNA sequencing (scRNA-seq) datasets from mouse and human samples, focusing on cross-species analysis through ortholog mapping and batch correction using Seurat and Harmony.

📜 Project Overview
This pipeline achieves the following objectives:

Data Preparation: Load and preprocess scRNA-seq datasets from human and mouse.
Ortholog Mapping: Map mouse genes to their human orthologs using a precomputed ortholog table.
Metadata Curation: Standardize metadata across datasets.
Integration: Merge datasets, normalize, scale, and perform PCA.
Batch Correction: Use Harmony to correct batch effects between species.
Visualization: Generate UMAP plots for cluster comparison across species.
