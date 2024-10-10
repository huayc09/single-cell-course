# Lesson 6: Advanced scRNA-seq Analytical Methods

# 1. SCENIC: Single-Cell Regulatory Network Inference and Clustering

library(Seurat)
library(SeuratExtend)

# Download the pre-computed SCENIC loom file
scenic_loom_path <- file.path(tempdir(), "pyscenic_integrated-output.loom")
download.file("https://zenodo.org/records/10944066/files/pbmc3k_small_pyscenic_integrated-output.loom", scenic_loom_path)

# Use the example PBMC dataset from SeuratExtend package
pbmc <- SeuratExtend::pbmc

# Import SCENIC results into our Seurat object
pbmc <- ImportPyscenicLoom(scenic_loom_path, seu = pbmc)

# View the AUCell matrix (TF activity scores)
tf_auc <- pbmc@misc$SCENIC$RegulonsAUC
head(tf_auc[, 1:5])

# View the list of TFs and their target genes
tf_gene_list <- pbmc@misc$SCENIC$Regulons
str(tf_gene_list, list.len = 10, max.level = 1)

# Visualizing SCENIC Results

# Identifying Top Activated TFs in Each Cluster
tf_zscore <- CalcStats(tf_auc, f = pbmc$cluster, order = "p", n = 4, t = TRUE)
Heatmap(tf_zscore, lab_fill = "zscore")

# Comparing TF Gene Expression Levels and Regulon Activity
DimPlot2(
  pbmc,
  features = c("ETS1", "ATF3", "tf_ETS1", "tf_ATF3"),
  cols = list("tf_ETS1" = "D", "tf_ATF3" = "D"),
  theme = NoAxes()
)

# Comparing Regulon Activity Between Cell Types
DefaultAssay(pbmc) <- "TF"
WaterfallPlot(
  pbmc,
  features = rownames(pbmc),
  ident.1 = "Mono CD14",
  ident.2 = "CD8 T cell",
  exp.transform = FALSE,
  top.n = 20
)

# 2. Trajectory Analysis: scVelo and Palantir

# Install miniconda and create conda environment (if not already done)
# reticulate::install_miniconda()
# SeuratExtend::create_condaenv_seuratextend()

library(Seurat)
library(SeuratExtend)

# Download the example Seurat Object and loom file
mye_small <- readRDS(url("https://zenodo.org/records/10944066/files/pbmc10k_mye_small_velocyto.rds", "rb"))
loom_path <- file.path(tempdir(), "pbmc10k_mye_small.loom")
download.file("https://zenodo.org/records/10944066/files/pbmc10k_mye_small.loom", loom_path)

# Prepare AnnData object for scVelo
adata_path <- file.path(tempdir(), "mye_small.h5ad")
scVelo.SeuratToAnndata(
  mye_small,
  filename = adata_path,
  velocyto.loompath = loom_path,
  prefix = "sample1_",
  postfix = "-1"
)

# Generate basic scVelo plot
scVelo.Plot(color = "cluster", save = "scvelo_basic.png", figsize = c(5,4))

# Generate a more customized scVelo plot
scVelo.Plot(
  style = "scatter",
  color = "cluster",
  groups = c("DC", "Mono CD14"),
  palette = color_pro(3, "light"),
  xlim = c(0, 5), ylim = c(0, 10),
  save = "scvelo_custom.png",
  figsize = c(5,4)
)

# Palantir Analysis

# Running Diffusion Map
mye_small <- Palantir.RunDM(mye_small)

# Visualize the multiscale space (ms) embedding
DimPlot2(mye_small, reduction = "ms", group.by = "cluster", label = TRUE)

# Calculating Pseudotime and Cell Fates
start_cell <- "sample1_GAGAGGTAGCAGTACG-1"
mye_small <- Palantir.Pseudotime(mye_small, start_cell = start_cell)
ps <- mye_small@misc$Palantir$Pseudotime

# Visualize pseudotime and cell fates
colnames(ps)[3:4] <- c("fate1", "fate2")
mye_small@meta.data[,colnames(ps)] <- ps
DimPlot2(mye_small, features = colnames(ps), reduction = "ms", 
         cols = list(Entropy = "D"), theme = NoAxes())

# Comparing Gene Expression Along Trajectories
GeneTrendCurve.Palantir(mye_small, features = c("CD14", "FCGR3A"), pseudotime.data = ps)

GeneTrendHeatmap.Palantir(
  mye_small, 
  features = c("CD14", VariableFeatures(mye_small)[1:10]), 
  pseudotime.data = ps, 
  magic = FALSE, 
  lineage = "fate1"
)

# 3. Cell-Cell Communication: CellChat and NicheNet

# CellChat Analysis
library(CellChat)
library(patchwork)

data.input <- GetAssayData(pbmc, assay = "RNA", slot = "data")
labels <- Idents(pbmc)
meta <- data.frame(group = labels, row.names = names(labels))

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "group")
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Visualizing CellChat Results
netVisual_circle(cellchat@net$count, vertex.weight = table(cellchat@idents), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_heatmap(cellchat)

# NicheNet Analysis
library(nichenetr)
library(tidyverse)

# Load pre-processed Seurat object
seuratObj <- readRDS(url("https://zenodo.org/record/3531889/files/seuratObj.rds"))
seuratObj <- UpdateSeuratObject(seuratObj)

organism <- "mouse"

if(organism == "human"){
  lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
  weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
} else if(organism == "mouse"){
  lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
  weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
}

# Define sender and receiver cells
sender_celltypes <- c("CD4 T", "Treg", "Mono", "NK", "B", "DC")
receiver_celltype <- "CD8 T"

# Run NicheNet analysis
nichenet_output <- nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj, 
  sender = sender_celltypes,
  receiver = receiver_celltype, 
  condition_colname = "aggregate",
  condition_oi = "LCMV",
  condition_reference = "SS",
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks
)

# View top predicted ligands
head(nichenet_output$top_ligands, 10)

# Visualize ligand activity and target genes
nichenet_output$ligand_activity_target_heatmap

# Visualize ligand expression in sender cells
nichenet_output$ligand_expression_dotplot

# 4. Label Transfer with singleCellNet

library(singleCellNet)
library(Seurat)
library(SingleCellExperiment)
library(tidyverse)

# Load your reference and query Seurat objects
reference_seurat <- pbmc
query_seurat <- readRDS("rds/pbmc_merge.rds")

# Convert Seurat objects to SingleCellExperiment objects
sce_reference <- as.SingleCellExperiment(reference_seurat)
sce_query <- as.SingleCellExperiment(query_seurat)

# Find common genes between the two datasets
commonGenes <- intersect(rownames(sce_reference), rownames(sce_query))

# Subset both datasets to only include common genes
sce_reference <- sce_reference[commonGenes, ]
sce_query <- sce_query[commonGenes, ]

# Prepare training data
stTrain <- data.frame(
  cell_name = rownames(colData(sce_reference)), 
  cell_type = colData(sce_reference)$cluster
)

expTrain <- assays(sce_reference)$counts

# Train the classifier
class_info <- scn_train(
  stTrain = stTrain, 
  expTrain = expTrain, 
  nTopGenes = 10, 
  nRand = 70, 
  nTrees = 1000, 
  nTopGenePairs = 25, 
  dLevel = "cell_type", 
  colName_samp = "cell_name"
)

# Classify cells in the query dataset
classRes_query <- scn_predict(
  cnProc = class_info[['cnProc']], 
  expDat = assays(sce_query)$counts, 
  nrand = 50
)

# Add classification results to the query Seurat object
query_seurat$Predicted_Labels <- rownames(classRes_query)[apply(classRes_query[,colnames(query_seurat)], 2, which.max)]

# Visualize the results
DimPlot2(query_seurat, group.by = "Predicted_Labels", label = TRUE)

# 5. MAGIC for Denoising and Smoothing Gene Expression

# Run MAGIC
mye_small <- Palantir.Magic(mye_small)

# Normalize the new "magic" assay
mye_small <- NormalizeData(mye_small)

# Visualize the effects of MAGIC
DimPlot2(mye_small, features = c("CD14", "magic_CD14", "FLT3", "magic_FLT3"), theme = NoAxes())

# 6. Copy Number Variation Analysis: CopyKat

# Install CopyKat if not already installed
# if (!require("copykat")) remotes::install_github("navinlabcode/copykat")

library(copykat)

# Assuming 'exp.rawdata' is your raw UMI count matrix
# copykat.test <- copykat(
#   rawmat=exp.rawdata, 
#   id.type="S", 
#   cell.line="no", 
#   ngene.chr=5, 
#   win.size=25, 
#   KS.cut=0.1, 
#   sam.name="test", 
#   distance="euclidean", 
#   n.cores=1)

# View prediction results
# head(copykat.test$prediction)

# View copy number matrix
# head(copykat.test$CNAmat[, 1:5])

# 7. TCR/BCR Analysis: scRepertoire

# Install scRepertoire if not already installed
# BiocManager::install("scRepertoire")

library(scRepertoire)

# Read the filtered_contig_annotations.csv files
# S1 <- read.csv("path/to/Sample1/outs/filtered_contig_annotations.csv")
# S2 <- read.csv("path/to/Sample2/outs/filtered_contig_annotations.csv")

# Create a list of the contig data
# contig_list <- list(S1, S2)

# Process the contig data
# combined <- combineTCR(contig_list, samples = c("sample1", "sample2"))
