# Lesson 2. Basic Single-Cell Analysis Workflow with Seurat

# This script demonstrates a basic single-cell RNA-seq analysis workflow using the Seurat package.

# 1. Load Required Libraries
# Make sure you have installed the Seurat package. If not, uncomment and run the following line:
# install.packages("Seurat")

library(Seurat)

# 2. Load the Data

# Set working directory to the course folder created in Lesson 1
# If you used a different location, replace this path with your chosen directory
setwd("~/Documents/single-cell-course")

# Download the course data (this may take a few minutes depending on your internet connection)
download.file(
  "https://zenodo.org/records/14245429/files/single-cell-course-data.zip",
  destfile = "course_data.zip"
)

# Extract the downloaded zip file
unzip("course_data.zip", exdir = "data")

# Using a 10x Genomics dataset of 3k PBMCs
raw_matrix <- Read10X("data/pbmc3k/filtered_gene_bc_matrices/")
raw_matrix[101:120,1:20]  # View a portion of the raw count matrix

# 3. Create a Seurat Object
pbmc <- CreateSeuratObject(
  raw_matrix,
  min.features = 200,  # Keep cells with at least 200 detected genes
  project = "pbmc_3k", # Name of the project
  min.cells = 1        # Keep genes detected in at least 1 cell
)

# Explore the Seurat object
View(pbmc)  # Opens a viewer to explore the object

# View a small portion of the raw count matrix
pbmc@assays$RNA$counts[1:20, 1:20]

# View metadata associated with cells
head(pbmc@meta.data)

# View cell names and gene names
head(colnames(pbmc), 10)  # View cell names
head(rownames(pbmc), 10)  # View gene names

# 4. Normalize the Data
pbmc <- NormalizeData(pbmc)

# View a portion of the normalized data
pbmc@assays$RNA$data[1:20, 1:10]

# 5. Identify Variable Features
pbmc <- FindVariableFeatures(pbmc)  # Identify genes that are highly variable across cells
head(VariableFeatures(pbmc), 20)  # View the variable features

# 6. Scale the Data and Perform PCA
pbmc <- ScaleData(pbmc)  # Scale and center the data
pbmc <- RunPCA(pbmc)

# View PCA results
head(Embeddings(pbmc, reduction = "pca"), 3)

# 7. Cluster the Cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)  # Construct a KNN graph
pbmc <- FindClusters(pbmc, resolution = 0.5)  # Find clusters of cells

# View cluster assignments
head(pbmc@meta.data)

# 8. Run UMAP for Visualization
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Visualize the clusters on the UMAP
DimPlot(pbmc, label = TRUE)

# 9. Annotate Cell Types
# Visualize expression of known marker genes
FeaturePlot(
  pbmc,
  features = c(
    "MS4A1", "GNLY", "CD3E",
    "CD8A", "CCR7", "CD14",
    "FCER1A", "FCGR3A", "PPBP"),
  order = TRUE
)

# Assign cell type labels
pbmc$cluster <- pbmc$seurat_clusters
head(pbmc$cluster)

new.cluster.ids <- c("CD4 T Naive", "CD4 T Memory", "Mono CD14", "B cell", "CD8 T cell", "Mono FCGR3A", "NK cell", "DC", "Platelet")
levels(pbmc$cluster) <- new.cluster.ids
head(pbmc$cluster)

# View the new annotation results
DimPlot(pbmc, label = TRUE, group.by = "cluster")

# Sort the levels alphabetically
pbmc$cluster <- factor(pbmc$cluster, levels = sort(levels(pbmc$cluster)))

# Visualize the final annotated UMAP
DimPlot(pbmc, label = TRUE, group.by = "cluster")

# Change the default idents
Idents(pbmc) <- 'cluster'
DimPlot(pbmc, label = TRUE)

# Save the Seurat object as an RDS file
dir.create("rds", showWarnings = FALSE)
saveRDS(pbmc, file = "rds/pbmc_annotated.rds")

# Install SeuratExtend for future use
# If you haven't installed remotes yet:
# install.packages("remotes")
remotes::install_github("huayc09/SeuratExtend")

# This concludes our basic single-cell analysis workflow with Seurat.
# In the next lessons, we'll explore more advanced techniques and visualizations using SeuratExtend.
