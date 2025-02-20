setwd("C:/Users/crtuser/OneDrive - TCDUD.onmicrosoft.com/Documents/PhD/Project/Data/WGCNA")


#---- scGNN
library(WGCNA)
library(hdWGCNA)
library(Seurat)
library(Matrix)
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggradar)
library(ggraph)
library(tidygraph)
library(fgsea)
library(dplyr)
library(magrittr)
library(STRINGdb)

# Download the file using R's built-in function
download.file("https://swaruplab.bio.uci.edu/public_data/hg38_mm10_orthologs_2021.txt", 
              destfile = "hg38_mm10_orthologs_2021.txt", 
              method = "auto")

# load in both files:
seurat_mouse <- readRDS("~/PhD/Project/Data/Hasty CiteSeq/GSE182233_IntegratedData.rds/GSE182233_IntegratedData.rds")
seurat_obj <- readRDS('C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\HumanWeightLossCohort\\WL_PATIENT_10X_Human\\D19-7153\\filtered_feature_bc_matrix\\200125_humanWL_cells.rds')

colnames(seurat_mouse@meta.data)[grep("orig.ident", colnames(seurat_mouse@meta.data), ignore.case = TRUE)] <- "Condition"


seurat_mouse <- subset(seurat_mouse, Condition == 'Lean', invert = T)
seurat_obj@meta.data$highlevel2 <- Idents(seurat_obj)

seurat_obj@meta.data$Condition <- ifelse(
  grepl("w1", seurat_obj@meta.data$orig.ident), "Obese",
  ifelse(grepl("w6", seurat_obj@meta.data$orig.ident), "WL", NA)
)

# Remove columns by index or name
seurat_obj@meta.data <- seurat_obj@meta.data[, !colnames(seurat_obj@meta.data) %in% c("celldex", "orig.ident")]

# Verify removal
head(seurat_obj@meta.data)


hg38_mm10_genes <- read.table("hg38_mm10_orthologs_2021.txt", sep='\t', header=TRUE)

colnames(hg38_mm10_genes) <-c('hg38_id', 'mm10_id', 'mm10_name', 'hg38_name')

# remove entries that don't have an ortholog
hg38_mm10_genes <- subset(hg38_mm10_genes, mm10_name != '' & hg38_name != '')

# show what the table looks like
head(hg38_mm10_genes)

mm10_genes <- unique(hg38_mm10_genes$mm10_name)
hg38_genes <- unique(hg38_mm10_genes$hg38_name)
hg38_mm10_genes <- hg38_mm10_genes[match(mm10_genes, hg38_mm10_genes$mm10_name),]

# get the mouse counts matrix, keep only genes with a human ortholog
X_mouse <- GetAssayData(seurat_mouse, layer ='counts')
X_mouse <- X_mouse[rownames(X_mouse) %in% hg38_mm10_genes$mm10_name,]

# rename mouse genes to human ortholog
ix <- match(rownames(X_mouse), hg38_mm10_genes$mm10_name)
converted_genes <- hg38_mm10_genes[ix,'hg38_name']
rownames(X_mouse) <- converted_genes
colnames(X_mouse) <- paste0(colnames(X_mouse), '_mouse')

# get the human counts matrix
X_human <- GetAssayData(seurat_obj, slot='counts')

# what genes are in common?
genes_common <- intersect(rownames(X_mouse), rownames(X_human))
X_human <- X_human[genes_common,]
colnames(X_human) <- paste0(colnames(X_human), '_human')

# make sure to only keep genes that are common in the mouse and human datasets
X_mouse <- X_mouse[genes_common,]

# set up metadata table for mouse
mouse_meta <- seurat_mouse@meta.data %>% dplyr::select(c(highlevel2, lowlevel2))
rownames(mouse_meta) <- colnames(X_mouse)

# set up metadata table for mouse
human_meta <- seurat_obj@meta.data %>% dplyr::select(c(Condition, highlevel2))
rownames(human_meta) <- colnames(X_human)

# make mouse seurat obj
seurat_m <- CreateSeuratObject(X_mouse, meta = mouse_meta)
seurat_h <- CreateSeuratObject(X_human, meta = human_meta)

# merge:
seurat_m$Species <- 'mouse'
seurat_h$Species <- 'human'

# Normalize and Scale individually BEFORE merging
seurat_m <- NormalizeData(seurat_m)
seurat_m <- FindVariableFeatures(seurat_m, nfeatures = 3000)
seurat_m <- ScaleData(seurat_m)

seurat_h <- NormalizeData(seurat_h)
seurat_h <- FindVariableFeatures(seurat_h, nfeatures = 3000)
seurat_h <- ScaleData(seurat_h)

seurat_merged <- merge(seurat_m, seurat_h)

# Identify variable features
common_features <- intersect(VariableFeatures(seurat_m), VariableFeatures(seurat_h))

# Ensure Harmony runs on common features
length(common_features)

seurat_object <- JoinLayers(seurat_merged)

rm(hg38_mm10_genes, human_meta, mouse_meta, seurat_h, seurat_m, seurat_mouse, seurat_obj, X_human, X_mouse)

seurat_merged <- NormalizeData(seurat_merged)
seurat_merged <- FindVariableFeatures(seurat_merged, nfeatures=3000)
seurat_merged <- ScaleData(seurat_merged)
seurat_merged <- RunPCA(seurat_merged)
# PCA Plot before Harmony
DimPlot(seurat_merged, reduction = "pca", group.by = "Species", label = TRUE)
seurat_merged <- RunHarmony(seurat_merged, group.by.vars = 'Species', theta = 0, lambda = 0.5)

seurat_merged <- RunUMAP(seurat_merged, reduction='harmony', dims=1:30, min.dist=0.3)

DimPlot(seurat_merged, group.by='highlevel2', split.by='Species', raster=FALSE, label=TRUE) + umap_theme() + NoLegend()