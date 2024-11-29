# Lesson 5: Core scRNA-seq Workflow Enhancements

# Load required libraries
library(Seurat)
library(SeuratExtend)

# Set working directory to the course folder created in Lesson 1
# If you used a different location, replace this path with your chosen directory
setwd("~/Documents/single-cell-course")

# 1. Merging Multiple Datasets

# Load 3' dataset
pbmc_3p <- Read10X("data/pbmc_ACDA/3p_ACDA/filtered_feature_bc_matrix/")
pbmc_3p <- CreateSeuratObject(pbmc_3p, project = "pbmc_3p", min.features = 200, min.cells = 1)

# Load 5' dataset
pbmc_5p <- Read10X("data/pbmc_ACDA/5p_ACDA/filtered_feature_bc_matrix/")
pbmc_5p <- CreateSeuratObject(pbmc_5p, project = "pbmc_5p", min.features = 200, min.cells = 1)

# Merge the datasets
pbmc_merge <- merge(pbmc_3p, pbmc_5p, add.cell.ids = c("pbmc_3p","pbmc_5p"))
pbmc_merge <- JoinLayers(pbmc_merge) # This step is required for Seurat v5
pbmc_merge

# 2. Detailed Quality Control (nGene, percent.mt)

# Examine nGene distribution
VlnPlot2(pbmc_merge, features = "nFeature_RNA", group.by = "orig.ident")

# Calculate percent.mt
pbmc_merge[["percent.mt"]] <- PercentageFeatureSet(pbmc_merge, pattern = "^MT-")
VlnPlot2(pbmc_merge, features = "percent.mt", group.by = "orig.ident")

# Apply quality control filters
print(paste("Before filtering:", ncol(pbmc_merge), "cells"))

pbmc_merge <- subset(
  pbmc_merge,
  subset = nFeature_RNA > 1000 &
    nFeature_RNA < 7500 &
    percent.mt < 20
)

print(paste("After filtering:", ncol(pbmc_merge), "cells"))

# 3. Doublet Removal

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("scDblFinder", quietly = TRUE))
  BiocManager::install("scDblFinder")
library(scDblFinder)

# Convert Seurat object to SingleCellExperiment
sce <- as.SingleCellExperiment(pbmc_merge)

# Run scDblFinder
sce <- scDblFinder(sce, samples = "orig.ident")

# Add scDblFinder results to our Seurat object
pbmc_merge$scDblFinder.class <- sce$scDblFinder.class

# View doublet predictions
table(pbmc_merge$scDblFinder.class)

# Remove predicted doublets
pbmc_merge <- subset(pbmc_merge, subset = scDblFinder.class == "singlet")

# Check the number of cells after doublet removal
print(paste("After doublet removal:", ncol(pbmc_merge), "cells"))

# 4. Data Integration

# Traditional Seurat workflow
pbmc_merge <- NormalizeData(pbmc_merge)
pbmc_merge <- FindVariableFeatures(pbmc_merge)
pbmc_merge <- ScaleData(pbmc_merge)
pbmc_merge <- RunPCA(pbmc_merge)
pbmc_merge <- RunUMAP(pbmc_merge, dims = 1:10)
pbmc_merge <- FindNeighbors(pbmc_merge, dims = 1:10)
pbmc_merge <- FindClusters(pbmc_merge, resolution = 0.5)

DimPlot2(pbmc_merge, features = c("seurat_clusters", "orig.ident", "CD3D", "CD14"), theme = NoAxes())

# Run Harmony
if (!require("harmony", quietly = TRUE))
  install.packages("harmony")
library(harmony)

pbmc_merge <- RunHarmony(pbmc_merge, group.by.vars = "orig.ident")

# Proceed with usual workflow using Harmony-corrected dimensions
pbmc_merge <- RunUMAP(pbmc_merge, dims = 1:10, reduction = "harmony")
pbmc_merge <- FindNeighbors(pbmc_merge, dims = 1:10, reduction = "harmony")
pbmc_merge <- FindClusters(pbmc_merge, resolution = 0.5)

DimPlot2(pbmc_merge, features = c("seurat_clusters", "orig.ident", "CD3D", "CD14"), theme = NoAxes())

# Compare PCA and Harmony outputs
head(Embeddings(pbmc_merge, reduction = "pca"), 3)
head(Embeddings(pbmc_merge, reduction = "harmony"), 3)

# Visualize PCA and Harmony dimensions
VlnPlot2(pbmc_merge, features = paste0("PC_", 1:9), group.by = "orig.ident")
VlnPlot2(pbmc_merge, features = paste0("harmony_", 1:9), group.by = "orig.ident")

# Save Seurat object
saveRDS(pbmc_merge, file = "rds/pbmc_merge.rds")

# 5. Cell-Cycle Scoring and Regression

# Load cell cycle genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Calculate cell cycle scores
pbmc_merge <- CellCycleScoring(pbmc_merge, s.features = s.genes, g2m.features = g2m.genes)

# View the results
head(pbmc_merge@meta.data)

# Example of regressing out cell cycle scores (not run)
# pbmc_merge <- ScaleData(pbmc_merge, vars.to.regress = c("S.Score", "G2M.Score"))

# Example of regressing out cell cycle scores and percent.mt (not run)
# pbmc_merge <- ScaleData(pbmc_merge, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))

# 6. Alternative Normalization Methods (SCTransform)

if (!requireNamespace("glmGamPoi", quietly = TRUE)) 
    BiocManager::install("glmGamPoi")
options(future.globals.maxSize = 1000 * 1024^2)  # Increase to 1000 MiB

# SCTransform workflow
pbmc_merge <- SCTransform(pbmc_merge)
pbmc_merge <- RunPCA(pbmc_merge)
pbmc_merge <- RunHarmony(pbmc_merge, group.by.vars = "orig.ident", theta = 5)
pbmc_merge <- RunUMAP(pbmc_merge, dims = 1:20, reduction = "harmony")
pbmc_merge <- FindNeighbors(pbmc_merge, dims = 1:20, reduction = "harmony")
pbmc_merge <- FindClusters(pbmc_merge, resolution = 0.5)

DimPlot2(pbmc_merge, features = c("seurat_clusters", "orig.ident", "CD3D", "CD14"), theme = NoAxes())

# 7. Complete Workflow Summary

# Traditional Workflow
{
  # Load and merge datasets
  pbmc_3p <- Read10X("data/pbmc_ACDA/3p_ACDA/filtered_feature_bc_matrix/")
  pbmc_3p <- CreateSeuratObject(pbmc_3p, project = "pbmc_3p", min.features = 500, min.cells = 1)
  
  pbmc_5p <- Read10X("data/pbmc_ACDA/5p_ACDA/filtered_feature_bc_matrix/")
  pbmc_5p <- CreateSeuratObject(pbmc_5p, project = "pbmc_5p", min.features = 500, min.cells = 1)
  
  pbmc_merge <- merge(pbmc_3p, pbmc_5p, add.cell.ids = c("pbmc_3p","pbmc_5p"))
  pbmc_merge <- JoinLayers(pbmc_merge)
  
  # Quality control
  pbmc_merge[["percent.mt"]] <- PercentageFeatureSet(pbmc_merge, pattern = "^MT-")
  pbmc_merge <- subset(
    pbmc_merge,
    subset = nFeature_RNA > 1000 &
      nFeature_RNA < 7500 &
      percent.mt < 20
  )
  
  # Doublet removal
  sce <- as.SingleCellExperiment(pbmc_merge)
  sce <- scDblFinder(sce, samples = "orig.ident")
  pbmc_merge$scDblFinder.class <- sce$scDblFinder.class
  pbmc_merge <- subset(pbmc_merge, subset = scDblFinder.class == "singlet")
  
  # Traditional normalization and scaling
  pbmc_merge <- NormalizeData(pbmc_merge)
  pbmc_merge <- FindVariableFeatures(pbmc_merge)
  
  # Calculate cell cycle scores
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  pbmc_merge <- CellCycleScoring(pbmc_merge, s.features = s.genes, g2m.features = g2m.genes)
  
  # Scale data, regressing out percent.mt
  pbmc_merge <- ScaleData(pbmc_merge, vars.to.regress = c("percent.mt"))
  
  # Dimension reduction and clustering
  pbmc_merge <- RunPCA(pbmc_merge)
  pbmc_merge <- RunHarmony(pbmc_merge, group.by.vars = "orig.ident")
  pbmc_merge <- RunUMAP(pbmc_merge, dims = 1:10, reduction = "harmony")
  pbmc_merge <- FindNeighbors(pbmc_merge, dims = 1:10, reduction = "harmony")
  pbmc_merge <- FindClusters(pbmc_merge, resolution = 0.5)
  
  # Visualization
  DimPlot2(pbmc_merge, features = c("seurat_clusters", "orig.ident", "CD3D", "CD14"), theme = NoAxes())
  
  return(pbmc_merge)
}

# Integrated Workflow using RunBasicSeurat
{
  # Load and merge datasets
  pbmc_3p <- Read10X("data/pbmc_ACDA/3p_ACDA/filtered_feature_bc_matrix/")
  pbmc_3p <- CreateSeuratObject(pbmc_3p, project = "pbmc_3p", min.features = 500, min.cells = 1)
  
  pbmc_5p <- Read10X("data/pbmc_ACDA/5p_ACDA/filtered_feature_bc_matrix/")
  pbmc_5p <- CreateSeuratObject(pbmc_5p, project = "pbmc_5p", min.features = 500, min.cells = 1)
  
  pbmc_merge <- merge(pbmc_3p, pbmc_5p, add.cell.ids = c("pbmc_3p","pbmc_5p"))
  pbmc_merge <- JoinLayers(pbmc_merge)
  
  # Doublet removal
  sce <- as.SingleCellExperiment(pbmc_merge)
  sce <- scDblFinder(sce, samples = "orig.ident")
  pbmc_merge$scDblFinder.class <- sce$scDblFinder.class
  pbmc_merge <- subset(pbmc_merge, subset = scDblFinder.class == "singlet")
  
  # Run the integrated workflow
  pbmc_merge <- RunBasicSeurat(
    pbmc_merge,
    spe = "human",
    nFeature_RNA.min = 1000,
    nFeature_RNA.max = 7500,
    percent.mt.max = 20,
    dims = 1:10,
    resolution = 0.5,
    vars.to.regress = "percent.mt",
    harmony.by = "orig.ident"
  )
  
  # Visualization
  DimPlot2(pbmc_merge, features = c("seurat_clusters", "orig.ident", "CD3D", "CD14"), theme = NoAxes())
  
  return(pbmc_merge)
}

# SCTransform Workflow
{
  # Load and merge datasets
  pbmc_3p <- Read10X("data/pbmc_ACDA/3p_ACDA/filtered_feature_bc_matrix/")
  pbmc_3p <- CreateSeuratObject(pbmc_3p, project = "pbmc_3p", min.features = 500, min.cells = 1)
  
  pbmc_5p <- Read10X("data/pbmc_ACDA/5p_ACDA/filtered_feature_bc_matrix/")
  pbmc_5p <- CreateSeuratObject(pbmc_5p, project = "pbmc_5p", min.features = 500, min.cells = 1)
  
  pbmc_merge <- merge(pbmc_3p, pbmc_5p, add.cell.ids = c("pbmc_3p","pbmc_5p"))
  pbmc_merge <- JoinLayers(pbmc_merge)
  
  # Quality control
  pbmc_merge[["percent.mt"]] <- PercentageFeatureSet(pbmc_merge, pattern = "^MT-")
  pbmc_merge <- subset(
    pbmc_merge,
    subset = nFeature_RNA > 1000 &
      nFeature_RNA < 7500 &
      percent.mt < 20
  )
  
  # Doublet removal
  sce <- as.SingleCellExperiment(pbmc_merge)
  sce <- scDblFinder(sce, samples = "orig.ident")
  pbmc_merge$scDblFinder.class <- sce$scDblFinder.class
  pbmc_merge <- subset(pbmc_merge, subset = scDblFinder.class == "singlet")
  
  # Calculate cell cycle scores
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  pbmc_merge <- CellCycleScoring(pbmc_merge, s.features = s.genes, g2m.features = g2m.genes)
  
  # SCTransform normalization
  pbmc_merge <- SCTransform(pbmc_merge, vars.to.regress = c("percent.mt"))
  
  # Dimension reduction and clustering
  pbmc_merge <- RunPCA(pbmc_merge)
  pbmc_merge <- RunHarmony(pbmc_merge, group.by.vars = "orig.ident", theta = 5)
  pbmc_merge <- RunUMAP(pbmc_merge, dims = 1:20, reduction = "harmony")
  pbmc_merge <- FindNeighbors(pbmc_merge, dims = 1:20, reduction = "harmony")
  pbmc_merge <- FindClusters(pbmc_merge, resolution = 0.5)
  
  # Visualization
  DimPlot2(pbmc_merge, features = c("seurat_clusters", "orig.ident", "CD3D", "CD14"), theme = NoAxes())
  
  return(pbmc_merge)
}
