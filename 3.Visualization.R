# Lesson 3. Advanced Visualization Techniques for scRNA-seq Data

# This script demonstrates advanced visualization techniques for single-cell RNA-seq data
# using Seurat and SeuratExtend packages.

# Load required libraries
library(Seurat)
library(SeuratExtend)
library(ggplot2)
library(cowplot)
library(viridis)
library(RColorBrewer)

# Set working directory to the course folder created in Lesson 1
# If you used a different location, replace this path with your chosen directory
setwd("~/Documents/single-cell-course")

# Load the PBMC object saved at the end of Lesson 2
pbmc <- readRDS("rds/pbmc_annotated.rds")

# 1. Create an Enhanced Dimensional Reduction Plot

# Basic usage
DimPlot2(pbmc)

# Visualizing different variables
DimPlot2(pbmc, features = c("cluster", "MS4A1", "CD14", "CD3D"))

# Splitting by variable
set.seed(42)
pbmc$sample <- sample(c("sample1", "sample1", "sample2"), size = ncol(pbmc), replace = TRUE)
table(pbmc$sample)  # Check the distribution
DimPlot2(pbmc, features = c("cluster", "CD14"), split.by = "sample", ncol = 1)

# Highlighting specific cells
b_cells <- colnames(pbmc)[pbmc$cluster == "B cell"]
DimPlot2(pbmc, cells.highlight = b_cells)

# Advanced customization
DimPlot2(
  pbmc,
  features = c("cluster", "sample", "CD14", "CD3D"),
  cols = list(
    "cluster" = "light",
    "CD14" = "D",
    "CD3D" = c("#EEEEEE", "black")
  ),
  theme = NoAxes())

# Adding labels and boxes
DimPlot2(pbmc, label = TRUE, box = TRUE, label.color = "black", repel = TRUE, theme = NoLegend())

# Simplifying labels with indices
DimPlot2(pbmc, index.title = "C", box = TRUE, label.color = "black")

# Adding UMAP Arrows
DimPlot2(
  pbmc, 
  features = c("cluster", "sample", "CD14", "CD3D"),
  theme = NoAxes()
) + theme_umap_arrows()

# 2. Simultaneous Display of Three Features on a Dimension Reduction Plot

# Using RYB system
FeaturePlot3(pbmc, color = "ryb", feature.1 = "CD3D", feature.2 = "CD14", feature.3 = "CD79A", pt.size = 0.5)

# Using RGB system
FeaturePlot3(pbmc, color = "rgb", feature.1 = "CD3D", feature.2 = "CD14", feature.3 = "CD79A", pt.size = 1)

# Batch visualization with FeaturePlot3.grid
FeaturePlot3.grid(pbmc, features = c("CD3D", "CD14", "CD79A", "FCGR3A", NA, "LYZ"), pt.size = 0.5)

# 3. Visualize Cluster Distribution in Samples

# Basic usage
ClusterDistrBar(origin = pbmc$sample, cluster = pbmc$cluster)

# Displaying absolute cell counts
ClusterDistrBar(origin = pbmc$sample, cluster = pbmc$cluster, percent = FALSE)

# Vertical bar plot
ClusterDistrBar(origin = pbmc$sample, cluster = pbmc$cluster, flip = FALSE, reverse_order = TRUE)

# Non-stacked bar plot
ClusterDistrBar(origin = pbmc$sample, cluster = pbmc$cluster, flip = FALSE, stack = FALSE)

# Exporting data matrix
data_matrix <- ClusterDistrBar(origin = pbmc$sample, cluster = pbmc$cluster, plot = FALSE)
print(data_matrix)

# 4. Create an Enhanced Violin Plot

# Basic violin plot with box plot and points
genes <- c("CD3D", "CD14", "CD79A")
VlnPlot2(pbmc, features = genes, ncol = 1)

# Hide points and outliers for a cleaner appearance
VlnPlot2(pbmc, features = genes, pt = FALSE, hide.outlier = TRUE, ncol = 1)

# Grouping by cluster and splitting each cluster by samples
VlnPlot2(pbmc, features = genes, group.by = "cluster", split.by = "sample")

# Filtering for certain subtypes and arranging plots in columns
cells <- colnames(pbmc)[pbmc$cluster %in% c("B cell", "Mono CD14", "CD8 T cell")]
VlnPlot2(pbmc, features = genes, group.by = "cluster", cells = cells)

# Adding statistical annotations using the Wilcoxon test
VlnPlot2(pbmc, features = genes, group.by = "cluster", cells = cells, 
         stat.method = "wilcox.test", hide.ns = TRUE)

# Restricting statistical comparisons and using t-test
VlnPlot2(pbmc, features = genes, group.by = "cluster", cells = cells, 
         stat.method = "t.test", comparisons = list(c(1,2), c(1,3)), hide.ns = FALSE)

# Adding Mean and Median Lines
lowExprGenes <- c("CCR7", "IL7R", "TCF7")
VlnPlot2(pbmc, 
         features = lowExprGenes, 
         show.mean = TRUE,      # Show mean and median lines
         mean_colors = c("red", "blue"),  # Colors for mean and median
         cols = "light",        # Light color scheme for better visibility
         ncol = 1)

# 5. Generate a Heatmap Plot

# Generate sample matrix
genes <- VariableFeatures(pbmc)
toplot <- CalcStats(pbmc, features = genes, method = "zscore", order = "p", n = 4)
head(toplot, 10)

# Basic heatmap
Heatmap(toplot, lab_fill = "zscore")

# Find markers for B cells
bcell.markers <- FindMarkers(pbmc, ident.1 = "B cell", logfc.threshold = 1, only.pos = TRUE)
head(bcell.markers)

# Find markers for all clusters
pbmc.markers <- FindAllMarkers(pbmc, logfc.threshold = 2, only.pos = TRUE)
head(pbmc.markers, 10)

# Modifying axis and labels
Heatmap(toplot, lab_fill = "zscore", plot.margin = margin(l = 30))

# Showing subset of gene names
genes <- VariableFeatures(pbmc)[1:500]
toplot2 <- CalcStats(pbmc, features = genes, method = "zscore", order = "p")
Heatmap(toplot2, lab_fill = "zscore", feature_text_subset = genes[1:20], expand_limits_x = c(-0.5, 12))

# Faceting the heatmap
gene_groups <- sample(c("group1", "group2", "group3"), nrow(toplot2), replace = TRUE)
Heatmap(toplot2, lab_fill = "zscore", facet_row = gene_groups) +
  theme(axis.text.y = element_blank())

# 6. Create Enhanced Dot Plots

# Basic Usage
genes <- VariableFeatures(pbmc)[1:10]
DotPlot2(pbmc, features = genes)

# Grouped Features
grouped_features <- list(
  "B_cell_markers" = c("MS4A1", "CD79A"),
  "T_cell_markers" = c("CD3D", "CD8A", "IL7R"),
  "Myeloid_markers" = c("CD14", "FCGR3A", "S100A8")
)
DotPlot2(pbmc, features = grouped_features)

# Split Visualization
DotPlot2(pbmc, features = genes, group.by = "cluster", split.by = "orig.ident", show_grid = FALSE)

# Alternative Split Visualization
# Using colors instead of borders for split groups
DotPlot2(pbmc, 
         features = genes, 
         group.by = "cluster", 
         split.by = "orig.ident", 
         split.by.method = "color", 
         show_grid = FALSE)

# Customizing Appearance
DotPlot2(pbmc, 
         features = grouped_features, 
         color_scheme = "BuRd", 
         border = FALSE,        # Remove dot borders
         show_grid = FALSE,     # Remove grid lines
         flip = TRUE)          # Flip coordinates

# 7. Generate a Waterfall Plot

# Basic waterfall plot
genes <- VariableFeatures(pbmc)[1:80]
WaterfallPlot(
  pbmc, group.by = "cluster", features = genes,
  ident.1 = "Mono CD14", ident.2 = "CD8 T cell", length = "logFC")

# Focusing on extremes
WaterfallPlot(
  pbmc, group.by = "cluster", features = genes,
  ident.1 = "Mono CD14", ident.2 = "CD8 T cell", length = "logFC",
  top.n = 20)

# 8. Create Volcano Plots

# Basic Usage
VolcanoPlot(pbmc, 
            ident.1 = "B cell",
            ident.2 = "CD8 T cell")

# Customizing Thresholds
VolcanoPlot(
  pbmc,
  ident.1 = "B cell",
  ident.2 = "CD8 T cell",
  x.threshold = 0.5,    # Log fold change threshold
  y.threshold = 2     # -log10(p-value) threshold
)

# 9. Tips: Choosing the Right Visualization
# See details in "3.Visualization.Rmd"

# 10. Explore Color Functions

# Discrete variables with "light" color scheme
plot_grid(
  DimPlot2(pbmc, cols = "light", theme = NoAxes() + NoLegend()),
  ClusterDistrBar(pbmc$sample, pbmc$cluster, cols = "light", flip = FALSE, border = "black", reverse_order = FALSE) +
    theme(axis.title.x = element_blank())
)

# Continuous variables with "A" color scheme
Heatmap(toplot, lab_fill = "zscore", color_scheme = "A")

# Continuous variables with "D" color scheme
WaterfallPlot(
  pbmc, group.by = "cluster", features = VariableFeatures(pbmc)[1:80],
  ident.1 = "Mono CD14", ident.2 = "CD8 T cell", length = "logFC",
  top.n = 20, color_theme = "D")

# Using RColorBrewer palettes
markers_genes <- c(
  "MS4A1", "GNLY", "CD3E",
  "CD8A", "CCR7", "CD14",
  "FCER1A", "FCGR3A", "PPBP")
DimPlot2(
  pbmc,
  features = markers_genes,
  cols = "OrRd",
  theme = NoAxes())

DimPlot2(
  pbmc,
  features = markers_genes,
  cols = "Spectral-rev",
  theme = NoAxes())

# 11. Mastering ggplot2 for Custom Visualizations

# Extract data from the Seurat object
umap_data <- FetchData(pbmc, vars = c("umap_1", "umap_2", "cluster", "CD3D", "sample"))
head(umap_data)

# Create the most basic plot
ggplot(umap_data, aes(x = umap_1, y = umap_2, color = CD3D)) +
  geom_point()

# Create an improved plot
ggplot(umap_data, aes(x = umap_1, y = umap_2, color = CD3D)) +
  geom_point(size = 0.5, alpha = 0.7) +
  scale_color_viridis() +
  theme_minimal() +
  labs(title = "UMAP colored by CD3D expression",
       x = "UMAP_1", y = "UMAP_2")

# Customizing the plot further
ggplot(umap_data, aes(x = umap_1, y = umap_2, color = CD3D)) +
  geom_point(size = 0.5, alpha = 0.7) +
  scale_color_viridis(option = "A") +
  facet_wrap(~sample, ncol = 3) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold")) +
  labs(title = "UMAP by cluster, colored by CD3D expression",
       x = "UMAP_1", y = "UMAP_2")

# Combining multiple plots
plot1 <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = cluster)) +
  geom_point(size = 0.5, alpha = 0.7) +
  scale_color_manual(values = color_pro(9,2)) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Clusters")

plot2 <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = CD3D)) +
  geom_point(size = 0.5, alpha = 0.7) +
  scale_color_viridis_c(option = "A") +
  theme_minimal() +
  labs(title = "CD3D Expression")

combined_plot <- plot_grid(plot1, plot2, labels = c("A", "B"), ncol = 2)

title <- ggdraw() + 
  draw_label("UMAP Visualization of PBMC Data", fontface = "bold", x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))

plot_grid(title, combined_plot, ncol = 1, rel_heights = c(0.1, 1))

# Saving the plot
dir.create("results", showWarnings = FALSE)
ggsave("results/pbmc_umap_plot.png", combined_plot, width = 8, height = 4.5, dpi = 300)

# This script provides advanced visualization techniques for single-cell RNA-seq data.
# In the next lesson, we'll explore more advanced analytical techniques to further your understanding of single-cell data.
