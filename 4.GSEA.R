# Lesson 4: Gene Set Enrichment Analysis and Utility Functions

# Load required libraries
library(Seurat)
library(SeuratExtend)
library(dplyr)
options(max.print = 15, spe = "human")

# Load the PBMC object we saved at the end of Lesson 2
pbmc <- readRDS("rds/pbmc_annotated.rds")

# Gene Ontology (GO) Database Analysis
pbmc <- GeneSetAnalysisGO(pbmc, parent = "immune_system_process", nCores = 4)
matr <- pbmc@misc$AUCell$GO$immune_system_process
matr <- RenameGO(matr, add_n_gene = FALSE)
head(matr, c(5,3))

# View commonly used GO categories
GeneSetAnalysisGO()

# Visualizing GO GSEA Results

# 1. Heatmap
toplot <- CalcStats(matr, f = pbmc$cluster, order = "p", n = 4)
Heatmap(toplot, lab_fill = "zscore")

# 2. Violin Plot
VlnPlot2(matr[1:3,], f = pbmc$cluster, ncol = 1)

# 3. Waterfall Plot
WaterfallPlot(matr, f = pbmc$cluster, ident.1 = "B cell", ident.2 = "CD8 T cell", top.n = 20)

# Reactome Database Analysis
pbmc <- GeneSetAnalysisReactome(pbmc, parent = "Immune System")
matr <- pbmc@misc$AUCell$Reactome$`Immune System`
matr <- RenameReactome(matr)
Heatmap(CalcStats(matr, f = pbmc$cluster, order = "p", n = 4), lab_fill = "zscore")

# View commonly used Reactome categories
GeneSetAnalysisReactome()

# Perform GSEA using Customized Gene Sets

# Using the Hallmark 50 Gene Set
str(hall50$human, list.len = 10, max.level = 1)

pbmc <- GeneSetAnalysis(pbmc, genesets = hall50$human)
matr <- pbmc@misc$AUCell$genesets
Heatmap(CalcStats(matr, f = pbmc$cluster), lab_fill = "zscore")

# Exploring Additional Gene Sets
names(Genesets_data$human$GSEA)
names(PanglaoDB_data$marker_list_human)

# Converting Pathway IDs to Names
RenameGO(c("GO:0002376","GO:0050896"), spe = "human")
RenameReactome(c("R-HSA-109582","R-HSA-112316"), spe = "human")

# Filtering Pathway Lists

# Filtering GO Pathways
terms <- FilterGOTerms(parent = "GO:0002376")
RenameGO(terms)

terms2 <- FilterGOTerms(term = terms, n.min = 10, n.max = 100)
RenameGO(terms2)

terms3 <- FilterGOTerms(term = terms, only.end.terms = TRUE)
RenameGO(terms3)

# Filtering Reactome Pathways
terms <- FilterReactomeTerms(parent = "R-HSA-168256")
RenameReactome(terms)

# Finding Pathways in Databases

# General Search
result <- SearchDatabase(c("CD3D", "metabolic"))
names(result)
glimpse(head(result))

# Species-Specific Search
result <- SearchDatabase("antigen pr", spe = "mouse")
names(result)
glimpse(head(result))

# Customizing Return Types
result <- SearchDatabase("metabolic", return = "ID")
result

result <- SearchDatabase("metabolic", return = "genelist")
str(result, list.len = 10, max.level = 1)

# Creating GSEA Plots
GSEAplot(
  pbmc, 
  ident.1 = "CD4 T Naive", 
  title = "INTERFERON_GAMMA_RESPONSE",
  geneset = hall50$human$HALLMARK_INTERFERON_GAMMA_RESPONSE
)

# Utility Functions for Gene Analysis

# Gene Naming Conversions
human_genes <- VariableFeatures(pbmc)[1:6]
print(human_genes)

HumanToMouseGenesymbol(human_genes)
HumanToMouseGenesymbol(human_genes, match = FALSE)

# Converting Gene Expression Matrices
human_matr <- GetAssayData(pbmc)[human_genes, 1:4]
print(human_matr)
HumanToMouseGenesymbol(human_matr)

# Other Gene Naming Conversion Functions
mouse_genes <- c("Cd14", "Cd3d", "Cd79a")
MouseToHumanGenesymbol(mouse_genes, match = FALSE)

human_genes <- c("PPBP", "LYZ", "S100A9", "IGLL5", "GNLY", "FTL")
GenesymbolToEnsembl(human_genes, spe = "human", keep.seq = TRUE)

GenesymbolToEnsembl(mouse_genes, spe = "mouse", keep.seq = TRUE)

EnsemblToGenesymbol(c("ENSG00000163736", "ENSG00000090382"), spe = "human", keep.seq = TRUE)

EnsemblToGenesymbol(c("ENSMUSG00000051439", "ENSMUSG00000032094"), spe = "mouse", keep.seq = TRUE)

# Assess Proportion of Positive Cells in Clusters
options(max.print = 50)

genes <- VariableFeatures(pbmc)[1:10]

# Default usage
proportions <- feature_percent(pbmc, feature = genes)
print(proportions)

library(RColorBrewer)
Heatmap(proportions, lab_fill = "Proportion of\npositive cell", color_scheme = brewer.pal(5, "OrRd"))

# Adjusting the Expression Threshold
proportions_above_2 <- feature_percent(pbmc, feature = genes, above = 2)
print(proportions_above_2)

Heatmap(proportions_above_2, lab_fill = "Proportion of\nexpression > 2", color_scheme = brewer.pal(5, "OrRd"))

# Targeting Specific Clusters
proportions_subset <- feature_percent(pbmc, feature = genes, ident = c("B cell", "CD8 T cell"))
print(proportions_subset)

# Proportions Relative to Total Cell Numbers
proportions_total <- feature_percent(pbmc, feature = genes, total = TRUE)
print(proportions_total)

# Logical Output for Expression
expressed_logical <- feature_percent(pbmc, feature = genes, if.expressed = TRUE, min.pct = 0.2)
print(expressed_logical)
