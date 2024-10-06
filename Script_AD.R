#Cypress HPC bash code
#singularity
#module load singularity/3.9.0
#module load R/4.1.1-intel
#source setup_cypress.sh 
#export SINGULARITY_IMAGES=/lustre/project/hdeng2/nwadiugwu/data/
#singularity shell -s /bin/bash $SINGULARITY_IMAGES/seurat_4.3.0.sif 


library(celldex)
library(ggplot2)
library(PISCES)
library(viper)
library(Seurat)
library(SingleR)
library(pheatmap)
library(circlize)
library(dplyr)
library(factoextra)
library(tibble)
library(cluster)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(Biobase)

source(file.path("/lustre/project/hdeng2/nwadiugwu/data/functions", "hands-on-1_heatmap-function.R"))
source(file.path("/lustre/project/hdeng2/nwadiugwu/data/functions", "tutorial-functions.R"))
source(file.path("/lustre/project/hdeng2/nwadiugwu/data/functions", "hands-on-1_clustering-functions.R"))


dirs <- list.dirs(path = '/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered', recursive = F, full.names = F)

for(x in dirs){
  name <- gsub('_filtered','', x)
  
  cts <- ReadMtx(mtx = paste0('/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/',x,'/snrna_counts.mtx'), feature.column = 1,
                 features = paste0('/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/',x,'/genes.csv'),
                 cells = paste0('/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/',x,'/barcodes_rna.csv'))
  
  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts))
}


seqad <- readRDS("/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/UCI.rds")




# QC & filtering -----------------------
ncol(seqad) 
mean(colSums(seqad)) 

View(seqad@meta.data)

# filtering
seqad <- subset(seqad, subset = nCount_RNA > 800 &
                  nFeature_RNA > 500)


### We will first normalize the data and generate the gene expression signature using the 
### SCTransform function
seqad_s <- SCTransform(seqad, conserve.memory = T, vars.to.regress = "percent.mt", verbose = FALSE)


#### Compute the Principal components to generate a 2D plot using the RunPCA function

# Scale the data
seqad_s <- ScaleData(seqad_s)

# Run PCA
seqad_s <- RunPCA(seqad_s, features = VariableFeatures(object = seqad_s))
seqad_s <- RunUMAP(seqad_s, dims = 1:30, verbose = FALSE)
seqad_s <- FindNeighbors(seqad_s, dims = 1:30, verbose = FALSE)
seqad_s <- FindClusters(seqad_s, verbose = FALSE)
DimPlot(seqad_s)

# Define a color palette with 24 distinct colors
my_color_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#FCCDE5")

# Use 'rep_my_color_palette' as the color palette for DimPlot
DimPlot(seqad_s, reduction = 'pca', dims = c(1, 2), cols = my_color_palette, pt.size = 1.5) + NoLegend()

#### Check the standard deviation explained by each principal component
ElbowPlot(seqad_s)

#### Plot the cells using two higher principal components (#15 and #16 as an example) 
DimPlot(seqad_s, reduction = 'pca', dims = c(15, 16), cols = my_color_palette, pt.size = 1.5) + NoLegend()

#### Perform UMAP dimensionality reduction using the RunUMAP function
seqad_s <- RunUMAP(seqad_s, dims = 1:15, verbose = FALSE, 
                   metric="correlation")
DimPlot(seqad_s, reduction = 'umap', cols = my_color_palette, pt.size = 1.5) + NoLegend()
##################################################################################################

###### Unsupevised Clustering Analysis ################################################
seqad_s <- FindNeighbors(seqad_s, dims = 1:15)
seqad_s <- FindClusters(seqad_s, verbose = TRUE, algorithm=1, 
                        resolution=seq(0.1,1,by=0.01))

top15_PCs<-as.data.frame(seqad_s$pca@cell.embeddings[,1:15])
euc_distance<-dist(top15_PCs, method="euclidean")

Louvain_SilScores<-Louvain_SilhouetteScore(seqad_s, euc_distance)

saveRDS(Louvain_SilScores, "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/Louvain_SilScores.rds")

#### Plot the average Silhouette Score values 
plot(seq(0.1,1,by=0.01), Louvain_SilScores, ylab="Slihouette Scores", xlab="Resolution Parameter", pch=18, 
     cex=1.5, main = "Louvain - Average Silhouette Scores", xaxt='n ', ylim = c(0.2,0.5))
axis(1, at = seq(0.1,1,by=0.01), las=2)

# Generate the silhouette plot and use the color_palette
fviz_silhouette(silhouette(as.numeric(seqad_s$SCT_snn_res.0.1), euc_distance),
                palette = my_color_palette)

seqad_s = saveRDS(seqad_s, "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/seqad_n2.rds")


seqad_s <- readRDS("/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/seqad_n2.rds")


# perform standard workflow steps to figure out if we see any batch effects --------
seqad_filtered <- NormalizeData(object = seqad_s)
seqad_filtered <- FindVariableFeatures(object = seqad_filtered)
seqad_filtered <- ScaleData(object = seqad_filtered)
seqad_filtered <- RunPCA(object = seqad_filtered)
ElbowPlot(seqad_filtered)
seqad_filtered <- FindNeighbors(object = seqad_filtered, dims = 1:20)
seqad_filtered <- FindClusters(object = seqad_filtered)
seqad_filtered <- RunUMAP(object = seqad_filtered, dims = 1:20)


# plot
p1 <- DimPlot(seqad_s, reduction = 'umap', group.by = 'Pathology.group')
p2 <- DimPlot(seqad_s, reduction = 'umap', group.by = 'Msex')
p3 <- DimPlot(seqad_s, reduction = 'umap', group.by = 'Subject')
p4 <- DimPlot(seqad_s, reduction = 'umap', group.by = 'Age',
              cols = my_color_palette)
#cols = c('red','green','blue'))

grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)


# perform integration to correct for batch effects ------
obj.list <- SplitObject(seqad_s, split.by = 'Pathology.group')
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}


# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)

# find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features)

# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)
seurat.integrated = saveRDS(seurat.integrated, "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/seurat_integrated2.rds")

seurat.integrated <- readRDS("/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/seurat_integrated2.rds")


seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:15)


p5 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Pathology.group')
p6 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Msex')
p7 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Subject')
p8 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Age',
              cols = my_color_palette)


grid.arrange(p5, p6, p7, p8, ncol = 1, nrow = 1)

grid.arrange(p1, p5, p2, p6, ncol = 2, nrow = 2)
grid.arrange(p3, p7, p4, p8, ncol = 2, nrow = 2)
grid.arrange(p1, p5, p2, p6, p3, p7, p4, p8, ncol = 3, nrow = 3)


# combine 'Pathology.group' and 'Msex' into a new variable called 'Group'
seurat.integrated$Group <- paste(seurat.integrated$Pathology.group, seurat.integrated$Subject, seurat.integrated$Msex, sep = "-")

# Create custom color palette with distinct colors
custom_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494")


p1.1 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Group') +
  scale_fill_manual(values = custom_colors)

# Use facet_grid to create a grid of plots with distinct colors
p1.1 + facet_grid(.~Group, scales = 'free')

grid.arrange(p1.1, ncol = 1, nrow = 1)


### Use ggplot to color the cells in the PCA and UMAP plots according
### to the clustering solution


color_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494")

# Using the custom_color_palette for the plots

UMAP_dataframe<-data.frame(seurat.integrated@reductions[["umap"]]@cell.embeddings)
PCA_dataframe<-data.frame(seurat.integrated@reductions[["pca"]]@cell.embeddings)
Clusters<-as.factor(as.numeric(seurat.integrated$SCT_snn_res.0.1))

# Step 3: Prepare grouping information (use either "Pathology.group" or "Patient" column)
grouping_column <- seurat.integrated$Pathology.group  
grouping_factor <- as.factor(grouping_column)
colors <- c("red", "blue", "green")

ggplot(PCA_dataframe,aes(x = PC_1, y = PC_2, color=Clusters)) + theme_bw()+
  geom_point(size=2)+ggtitle("PCA Plot - Louvain Clusters")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  scale_color_manual(values= color_palette)

ggplot(UMAP_dataframe,aes(x = UMAP_1, y = UMAP_2, color=Clusters)) + theme_bw()+
  geom_point(size=2)+ggtitle("UMAP Plot - Louvain Clusters")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  scale_color_manual(values= color_palette)

blueprint.encode <- BlueprintEncodeData()
SingleR_Results <- SingleR(test = seurat.integrated$SCT@counts, ref = blueprint.encode, 
                           labels = blueprint.encode$label.main, assay.type.test = 2)


SingleR_Labels <- SingleR_Results$labels
SingleR_Labels[which(SingleR_Labels %in% names(which(table(SingleR_Labels)<10)))] <- "Others"

ggplot(seurat.integrated,aes(x = PC_1, y = PC_2, color=SingleR_Labels)) + theme_bw()+
  geom_point(size=2)+ggtitle("PCA Plot - SingleR Results")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  scale_color_manual(values=my_color_palette) 

ggplot(UMAP_dataframe,aes(x = UMAP_1, y = UMAP_2, color=SingleR_Labels)) + theme_bw()+
  geom_point(size=2)+ggtitle("UMAP Plot - SingleR Results")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  scale_color_manual(values=color_palette) 



######## Differential Gene Expression Analysis #################

#### Use the FindAllMarkers function of the Seurat pipeline to 
#### identify the most differnetially expressed genes per cluster through a Wilcoxon Test
Idents(seurat.integrated)<-"SCT_snn_res.0.1"
AD_markers <- FindAllMarkers(seurat.integrated, test.use = "wilcox", 
                             only.pos = TRUE, min.pct = 0.25)

#### Select the top 10 differentially expressed genes per cluster 
top10 <- AD_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

#### Plot the expression of the top10 genes in a heatmap.
### In the script Heatmap_Function.R there is the implementation of a function
### to generate and save in a .pdf file the figure of the heatmap

### Prepare the expression matrix to plot the genes
genes<-top10$gene
expression_to_plot=seurat.integrated[["SCT"]]@data[genes,]
Clusters<-as.factor(as.numeric(seurat.integrated$SCT_snn_res.0.1))
colors<- color_palette #c("red", "blue", "green", "purple", "orange", "cyan", "grey")

gExprHeatmap(genes, expression_to_plot, Clusters, colors)

write.csv(top10, file = "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/top10_gene_msex.csv", row.names = FALSE)


# Step 1: Extract the gene names from top10
genes <- top10$gene

# Step 2: Create the expression_to_plot matrix
expression_to_plot <- seurat.integrated[["SCT"]]@data[genes,]

# Step 3: Prepare grouping information (use either "Pathology.group" or "Patient" column)
#grouping_column <- seurat.integrated$Msex  # or use seurat.integrated$Patient if that's the grouping column
grouping_column <- paste(seurat.integrated$Msex, seurat.integrated$Age, seurat.integrated$Pathology.group, sep = "-")

# Step 4: Convert the grouping_column to a factor (if it's not already)
grouping_factor <- as.factor(grouping_column)

# Step 5: Create a vector of three distinct colors
#colors <- c("red", "blue", "green")
colors <- c("red", "blue", "green", "orange", "purple", "cyan")  # Add more colors if needed

# Step 6: Create the heatmap using gExprHeatmap with the grouping factor and the three colors
gExprHeatmap(genes, expression_to_plot, grouping_factor, colors)









# Part 2
#The part *Part A* will explore how to preprocess gene expression data to generate gene regulatory networks with ARACNe. In another part, *Part B*, we will run an analysis at the protein activity level.

### Import libraries and functions

rm(list = ls()) # cleans environment
cat("\014") # clean console

# install cran packages
#install.packages("abind", "BiocManager", "circlize", "cluster", "devtools",   "ggplot2", "ggpubr", "ggrepel", "grDevices", "Matrix",                  "RColorBrewer", "RSpectra", "Seurat", "uwot")

# install bioconductor packages
#BiocManager::install("biomaRt")
#BiocManager::install("ComplexHeatmap")
# install PISCES
#devtools::install_github("califano-lab/PISCES")
.libPaths()

suppressPackageStartupMessages({
  library(dplyr)
  library(Rcpp)
  library(harmony)
  library(Seurat)
  library(patchwork)
  library(qs)
  library(celldex)
  library(ggplot2)
  library(PISCES)
  library(viper)
  library(Seurat)
  library(SingleR)
  library(pheatmap)
  library(circlize)
  library(dplyr)
  library(factoextra)
  library(tibble)
  library(cluster)
  library(tidyverse)
  library(gridExtra)
  library(ggpubr)
  library(Biobase)
  ## Load pre-implemented functions
})


# Part A - Generating Gene Regulatory Networks

#**Goal**: to identify coarse cell groupings for network generation. 

#r Preprocessing at gene expression

# data
seqad_s <- readRDS("/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/seqad_n2.rds")


# Scale the data
seqad_s <- SCTransform(seqad_s, conserve.memory = TRUE, verbose = FALSE)

# Run PCA
seqad_s <- RunPCA(seqad_s, features = VariableFeatures(object = seqad_s))

# Extract the current row names from the meta.data DataFrame
original_row_names <- rownames(seqad_s@meta.data)

# Modify the row names by removing everything after the first comma in each row name
# This will split each row name at the first comma and retain only the first part
modified_row_names <- sapply(original_row_names, function(name) {
  split_name <- unlist(strsplit(name, ",", fixed = TRUE))
  return(split_name[1])
})

# Update the row names in the meta.data with the modified values
rownames(seqad_s@meta.data) <- modified_row_names

# Display the first few rows to verify the change
head(seqad_s@meta.data)

seqad_s <- ScaleData(seqad_s)

seqad_s <- RunPCA(seqad_s, verbose = FALSE)
seqad_s <- RunUMAP(seqad_s, dims = 1:15, verbose = FALSE)
seqad_s <- FindNeighbors(seqad_s, dims = 1:15, verbose = FALSE)
seqad_s <- FindClusters(seqad_s, verbose = FALSE)

DimPlot(seqad_s)
ncol(seqad_s) # 61472
mean(colSums(seqad_s)) #9645.539


#r Semi-supervised cell-type inference by SingleR}

blueprint.encode <- celldex::BlueprintEncodeData() # load reference

bp.singleR <- SingleR(test = seqad_s@assays$SCT@data, ref = blueprint.encode, labels = blueprint.encode$label.main, assay.type.test = 2)

# Assuming 'bp.singleR$labels' is a named vector or list
# Find unique labels in the 'bp.singleR$labels' object
unique_labels <- unique(bp.singleR$labels)

# Print the unique labels
print(unique_labels)

# Count the occurrences of each unique label
label_counts <- table(bp.singleR$labels)

# Print the total number of each unique label
print(label_counts)


#[Step 3:]{.underline} we group cells together based on their lineage (coarse grain assignment)

lymphoid.cell.types <- c('B-cells', 'CD4+ T-cells', 'CD8+ T-cells', 'NK cells')
neuronal.cell.types <- c('Neurons')
myeloid.cell.types <- c('Monocytes', 'Macrophages')
lymphatic.cell.types <- c('Endothelial cells')
connective.tissue.cell.types <- c('Chondrocytes', 'Adipocytes', 'Fibroblasts')
other.cell.types <- c('Others')

# Select the column with label "Barcode"
barcode_column <- seqad_s@meta.data[,"Barcode"]
#barcode_column <- seqad_s@meta.data[,"sample"]

# View the selected column
head(barcode_column)

lymphoid.cells <- barcode_column[which(bp.singleR$labels %in% lymphoid.cell.types)]
myeloid.cells <- barcode_column[which(bp.singleR$labels %in% myeloid.cell.types)]
neuronal.cell <- barcode_column[which(bp.singleR$labels %in% neuronal.cell.types)]
connective.tissue.cell <- barcode_column[which(bp.singleR$labels %in% connective.tissue.cell.types)]
lymphatic.cell <- barcode_column[which(bp.singleR$labels %in% lymphatic.cell.types)]
other.cell <- barcode_column[which(bp.singleR$labels %in% other.cell.types)]

# we prepare matrices for network generation. The goal is to preprocess matrices to be used as input to ARACNe-AP or ARACNe3.
#r Metacells generation

print(colnames(seqad_s@meta.data))

# construct a distance matrix
mta.dist <- dist(seqad_s@reductions$pca@cell.embeddings[,1:15])
#write.table(mta.dist, file = "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/mta.dist.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
#mta.dist = read.table("/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/mta.dist.txt")

# Get the filtered expression matrix using the GetAssayData function
filtered_expression_matrix <- Seurat::GetAssayData(object = seqad_s, assay = "RNA")

# Alternatively, we can access the filtered expression matrix directly from the counts slot:
filtered_expression_matrix <- seqad_s@assays$RNA@counts

colnames(filtered_expression_matrix) <- barcode_column

# Convert filtered_expression_matrix to a data frame
filtered_expression_df <- as.data.frame(filtered_expression_matrix)

# Save the new data frame to a file
#write.table(filtered_expression_df, file = "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/filtered_expression_matrix_as_dataframe.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
#filtered_expression_df <- read.table("/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/filtered_expression_matrix_as_dataframe.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

# Data manipulation
# Extract the first item before the comma from the existing row names
new_row_names <- sub(",.*", "", rownames(seqad_s@meta.data))

# Assign the modified row names to the rownames of the Seurat object
rownames(seqad_s@meta.data) <- new_row_names

# Save the modified Seurat object to a new file
#saveRDS(seqad_s, file = "/work/nwadiugw/seqad/seqad_filtered/seqad_n2_modified.rds")

# Print the modified Seurat object's metadata
print(seqad_s@meta.data)

#head(seqad_n$seurat_clusters)

# Example adjustment if clust.samps calculation is off
####}

## generate meta cells
meta.cells <- make_metacells(filtered_expression_df, mta.dist, seqad_s$seurat_clusters,
                             num.neighbors = 5, min.samps = 1000)       

## generate meta cells
meta.cells <- make_metacells(filtered_expression_df, mta.dist, neuronal_seurat_subset$seurat_clusters,
                             num.neighbors = 5, min.samps = 1000)
# Save the results (optional)
saveRDS(meta.cells, file = "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/new_meta_cells.rds")

#we can subset cell populations and create their metacells
#neuronal_seurat_subset = readRDS("/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/neuronal_seurat_subset.rds")
#myeloid_seurat_subset =readRDS("/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/myeloid_seurat_subset.rds")      
#meta.cells = readRDS("/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/meta_cells2.RDS")


meta.cells[[1]] %>% head()


# Check if meta.cells is empty or has issues
head(print(meta.cells))


# Print cluster labels and first few rows of data
print(head(seqad_s$seurat_clusters))
print(head(meta.cells))


# Convert meta.cells to data frame
meta_cells_df <- as.data.frame(meta.cells[[1]])

# Add rownames as a new column
meta_cells_df$gene <- rownames(meta_cells_df)

# Reorder columns to have CellID as the first column
meta_cells_df <- meta_cells_df %>%
  select(gene, everything())

# Remove rownames from the data frame
rownames(meta_cells_df) <- NULL


# Extract the current row names from the meta.data DataFrame
original_row_names <- colnames(meta_cells_df)

# Modify the row names by removing everything after the first comma in each row name
# This will split each row name at the first comma and retain only the first part
modified_row_names <- sapply(original_row_names, function(name) {
  split_name <- unlist(strsplit(name, ",", fixed = TRUE))
  return(split_name[1])
})

# Update the row names in the meta.data with the modified values
colnames(meta_cells_df) <- modified_row_names

# Display the first few rows to verify the change
head(colnames(meta_cells_df))

# Write the data frame to a tab-delimited text file
write.table(meta_cells_df, file = "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/Specific_Metacells/meta_cells_new1.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



#creating metacells for specific cell population

# Ensure the Seurat package is loaded
library(Seurat)

# Assuming neuron_cells contains the barcodes of neuron cells of interest
# And seqad_s is the full Seurat object you want to subset

# Convert the rownames (barcodes) from neuron_cells into a vector
neuron_barcodes <- neuronal.cell

# Use the Subset function to keep only the cells with matching barcodes in seqad_s
neuronal_seurat_subset <- subset(seqad_s, cells = neuron_barcodes)

# Check the result to ensure the subset operation worked as expected
head(print(neuronal_seurat_subset))

# Save the neuron_seurat_subset Seurat object to an RDS file
saveRDS(neuronal_seurat_subset, file = "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/neuronal_seurat_subset.rds")

meta.cells <- make_metacells(filtered_expression_df, mta.dist, neuronal_seurat_subset$seurat_clusters,
                             num.neighbors = 5, min.samps = 1000)

saveRDS(meta.cells, file = "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/meta_cells_neuronal2.RDS")






### ARACNe

#Metacells can be used as inputs for ARACNe. Please notice that running ARACNe-AP or ARACNe3 for scRNA-seq data requires an HPC to generate the networks. The script used for HPC can be found in the califano-lab GitHub.

# Part B - Analyzing scRNA-seq data at the Protein Activity Level

#[Step 1:] We generate a gene expression signature to be used as input to VIPER.

#Seurat SCTransform is used to generate the signature. Alternatively, you can use the PISCES' built-in method (Optional) as: `ges.mat <- internal_ges(pbmc.filt, norm.method = 'pflpf', est.method = 'map')`

#r Gene expression signature for VIPER

ges.mat <- as.matrix(seqad_s@assays$SCT@scale.data)


#[Step 2:]{.underline} We load and explore ARACNe-inferred gene regulatory networks.

#Lineage specific networks (myeloid and neurons, respectively) and the metaVIPER version of the VIPER algorithm will be used to infer protein activity. 
#metaVIPER is particularly useful when lineage relationships are not evident. The implementation of the ARACNe algorithm used for the generation of these networks is ARACNe-AP.

#```{r Loading and inspecting ARACNe networks}

expression_matrix <- Seurat::GetAssayData(object = seqad_s, assay = "RNA")
expression_file <- "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/expression_matrix.txt"
#write.table(expression_matrix, file = expression_file, sep = "\t", quote = FALSE)

# Load ARACNe output and expression matrix
aracne_output_file <- "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/Networks/ADnvtfs_network.txt"
expression_matrix <- read.table(expression_file, header = TRUE, row.names = 1, sep = "\t")

# Convert the expression matrix to ExpressionSet
eset <- ExpressionSet(assayData = as.matrix(expression_matrix))  #
regulon <- aracne2regulon(aracne_output_file, eset)
neuron <- pruneRegulon(regulon, cutoff=50, adaptive=FALSE)

# Include regulons to a list for metaViper
net.list <- list('neurons' = neuron)

PA <- viper(ges.mat, net.list, method = "none")

# Include VIPER analysis into the Seurat object
seqad_s[["VIPER"]] <- CreateAssayObject(counts=PA) # it will be stored as "sparse"
DefaultAssay(seqad_s) <- "VIPER"

## compute distance matrix 
viper.pca <- fast_pca(as.matrix(seqad_s[["VIPER"]]@counts), num.pcs = 18)   
viper.dist <- dist(viper.pca$x)
viper.dist_matrix <- as.matrix(viper.dist) 

viper.clust <- louvain_k(viper.dist, kmin = 5, kmax = 200, kstep=1)


# Split the data into subsets (adjust as needed). Done due to large vector size.
num_subsets <- 2
subset_size <- ceiling(nrow(viper.dist_matrix) / num_subsets)
subset_indices <- split(1:nrow(viper.dist_matrix), rep(1:num_subsets, each = subset_size, length.out = nrow(viper.dist_matrix)))

# Initialize an empty list to store clustering results
clustering_results <- list()

# Iterate over subsets and perform clustering
for (subset_idx in seq_along(subset_indices)) {
  subset_indices_current <- subset_indices[[subset_idx]]
  viper.dist_subset <- viper.dist_matrix[subset_indices_current, subset_indices_current]
  
  # Perform clustering using louvain_k or other methods
  viper.clust_subset <- louvain_k(viper.dist_subset, kmin = 5, kmax = 200, kstep = 1)
  
  # Store the results in the list
  clustering_results[[subset_idx]] <- viper.clust_subset
}
# Combine the clustering results from all subsets
combined_clusters <- unlist(lapply(clustering_results, function(result) result$opt.clust))

# Save the combined clustering results
saveRDS(combined_clusters, "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/Networks/combined_clusters_ADcotfs.rds")

# Load the combined clustering results
combined_clusters <- readRDS("/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/Networks/combined_clusters_Nvcotfs2_network.txt")  




top15_PCs<-as.data.frame(seqad_s$pca@cell.embeddings[,1:11])
euc_distance<-dist(top15_PCs, method="euclidean")

Louvain_SilScores<-Louvain_SilhouetteScore(seqad_s, euc_distance)

#### Plot the average Silhouette Score values 
plot(seq(0.1,1,by=0.01), Louvain_SilScores, ylab="Slihouette Scores", xlab="Resolution Parameter", pch=18, 
     cex=1.5, main = "Louvain - Average Silhouette Scores", xaxt='n ', ylim = c(0.2,0.5))
axis(1, at = seq(0.1,1,by=0.01), las=2)

my_color_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#FCCDE5")

my_color_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5")

# Generate the silhouette plot and use the color_palette
fviz_silhouette(silhouette(as.numeric(seqad_s$SCT_snn_res.0.8), euc_distance),
                palette = my_color_palette)

saveRDS(fviz_silhouette, "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/Networks/_New_Nvcotfs2_fviz_silhouette.rds")



# Master Regulator Analysis

# Add optimal clustering solution to the Seurat object
seqad_s <- AddMetaData(object = seqad_s, metadata = combined_clusters, col.name = 'viper.clust')

Idents(seqad_s) <- "viper.clust" # set identity for the Seurat object

# Generate UMAP embedding
seqad_s <- RunUMAP(seqad_s, features=rownames(seqad_s), metric="correlation", verbose=FALSE)
viper.umap <- seqad_s[["umap"]]@cell.embeddings

plot.df <- data.frame('UMAP1' = viper.umap[,1], 'UMAP2' = viper.umap[,2],
                      'Cluster' = as.factor(combined_clusters)) %>%
  rownames_to_column(., var="cell_id")

umap_plot <- ggplot(plot.df, aes(UMAP1, UMAP2)) +
  geom_point(aes(color = Cluster)) +
  ggtitle('Protein Activity Clustering - UMAP') +
  theme_bw()

# Save UMAP plot as PDF
ggsave("/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/umap_plot.pdf", umap_plot)


## Compare clusters at protein activity with SingleR annotations
SingleR_Labels <- bp.singleR %>% 
  as.data.frame() %>% 
  dplyr::select(labels) %>%
  rownames_to_column(., var = "cell_id")

SingleR_Labels$labels[which(SingleR_Labels$labels %in% names(which(table(SingleR_Labels$labels) < 10)))] <- "Others"  

## plot.df for plotting, including SingleR labels
plot.df <- left_join(x = plot.df, y = SingleR_Labels, by = "cell_id")

umap_singleR_plot <- ggplot(plot.df, aes(x = UMAP1, y = UMAP2, color = labels)) +
  theme_bw() +
  geom_point(size = 2) +
  ggtitle("UMAP Protein Activity - SingleR Results")

# Save UMAP plot with SingleR annotations as PDF
ggsave("/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/umap_singleR.pdf", umap_singleR_plot)


# Master Regulator Analysis
viper.markers <- StoufferClusters(as.matrix(seqad_s[["VIPER"]]@counts), clusters = as.factor(seqad_s$Pathology.group)) # for other grouping - viper.markers1 <- StoufferClusters(as.matrix(seqad_s[["VIPER"]]@counts), clusters = as.factor(seqad_s$viper.clust))

write.table(viper.markers, file = "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/Anderson/Results/viper.markers_pathology_New_Nvcotfs2.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


viper.markers <- StoufferClusters(as.matrix(seqad_s[["VIPER"]]@counts), clusters = as.factor(seqad_s$Diagnosis)) # for other grouping - viper.markers1 <- StoufferClusters(as.matrix(seqad_s[["VIPER"]]@counts), clusters = as.factor(seqad_s$viper.clust))
viper.markers

#Msex
viper.markers <- StoufferClusters(as.matrix(seqad_s[["VIPER"]]@counts), clusters = as.factor(seqad_s$Sex)) # for other grouping - viper.markers1 <- StoufferClusters(as.matrix(seqad_s[["VIPER"]]@counts), clusters = as.factor(seqad_s$viper.clust))
viper.markers

write.table(viper.markers, file = "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/Anderson/Results/viper.markers_sex_myeloid_ba46_cotfs.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

## Select the top regulators 
regulators <- Reduce(union, apply(viper.markers, 2, function(x) { 
  x <- sort(x, decreasing = TRUE); 
  names(c(head(x, 5))) 
}))



# Assuming `viper.markers` and `seqad_s$Age` are correctly defined and loaded

# Select the top regulators
regulators <- Reduce(union, apply(viper.markers, 2, function(x) { 
  x <- sort(x, decreasing = TRUE)
  names(c(head(x, 5))) 
}))

# Create empty lists for upregulated and downregulated regulators
upregulated_regulators <- list()
downregulated_regulators <- list()

# Ensure the cluster names in `seqad_s$Diagnosis` match the column names in `viper.markers`
cluster_names <- intersect(unique(as.character(seqad_s$braak)), colnames(viper.markers))

# Loop through each cluster in viper.markers
for (cluster_name in cluster_names) {  
  # Get the scores for the current cluster
  x <- viper.markers[, cluster_name]  # Access cluster by name
  
  # Sort the scores in descending order
  sorted_x <- sort(x, decreasing = TRUE)
  
  # Ensure there are enough elements to select
  num_elements <- length(sorted_x)
  top_n <- min(5, num_elements)
  
  # Extract the top 5 upregulated regulators (or fewer if there are not enough elements)
  top_upregulated <- names(sorted_x[1:top_n])
  
  # Extract the top 5 downregulated regulators (or fewer if there are not enough elements)
  # Ensure the downregulated ones do not overlap with the upregulated ones
  if (num_elements > top_n) {
    top_downregulated <- names(sorted_x[(num_elements - top_n + 1):num_elements])
    top_downregulated <- setdiff(top_downregulated, top_upregulated)
  } else {
    top_downregulated <- character(0)  # No downregulated if there are not enough elements
  }
  
  # Store the results in the respective lists
  upregulated_regulators[[cluster_name]] <- top_upregulated
  downregulated_regulators[[cluster_name]] <- top_downregulated
}



# List of upregulated regulators for each cluster
print(upregulated_regulators)

# List of downregulated regulators for each cluster
print(downregulated_regulators)



# Save upregulated regulators to a text file
write.table(regulators, file = "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/Anderson/Results/regulators_DE_Sex_myeloid_ba46_tfs.tfs", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# Save upregulated regulators to a text file
write.table(upregulated_regulators, file = "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/Anderson/Results/upregulated_regulators_DE_Sex_myeloid_ba46_tfs.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# Save downregulated regulators to a text file
write.table(downregulated_regulators, file = "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/Anderson/Results/downregulated_regulators_DE_Sex_myeloid_ba46_tfs.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)




#differentially expressed PA
viper.markers <- StoufferClusters(as.matrix(seqad_s[["VIPER"]]@counts), clusters = as.factor(seqad_s$Pathology.group))

write.table(viper.markers, file = "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/Results/viper.markers_pathology_Nvcotfs2.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

#saveRDS(viper.markers, "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/Results/viper.markers_pathology_ADnv_conetwork.rds")

# Identify the top differentially expressed regulators (upregulated or downregulated)
regulators <- Reduce(union, apply(viper.markers, 2, function(x) { 
  x <- sort(abs(x), decreasing = TRUE); # Sort by absolute association scores
  names(c(head(x, 5))) 
}))

# Create an empty data frame to store the results
changing_regulators_df <- data.frame(
  Regulators = character(0),
  `AD` = numeric(0),
  `Control` = numeric(0)
)

# Loop through each regulator (row) in viper.markers
for (regulator_name in rownames(viper.markers)) {
  
  # Get the association scores for the current regulator across all clusters
  x <- as.vector(viper.markers[regulator_name, ])
  
  # Check if the regulator changes from positive to negative or vice versa across all clusters
  if (sum(x > 0) > 0 && sum(x < 0) > 0) {
    # Append the regulator name to the 'Regulators' column
    regulators <- regulator_name
    
    # Append the association scores to the respective columns
    AD <- x[1]
    Control <- x[2]
    
    # Combine all the data into the data frame
    changing_regulators_df <- rbind(changing_regulators_df, data.frame(
      Regulators = regulators,
      `AD` = AD,
      `Control` = Control
    ))
  }
}

# Print the resulting data frame
print(changing_regulators_df)
#saveRDS(changing_regulators_df$Regulators, "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/Results/changing_regulators_pathology_New_ntfs11.rds")



#differentially expressed PA for amyloid group
viper.markers <- StoufferClusters(as.matrix(seqad_s[["VIPER"]]@counts), clusters = as.factor(seqad_s$Sex))
write.table(viper.markers, file = "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/Results/viper.markers_msex_Nvcotfs2.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# Identify the top differentially expressed regulators (upregulated or downregulated)
regulators <- Reduce(union, apply(viper.markers, 2, function(x) { 
  x <- sort(abs(x), decreasing = TRUE); # Sort by absolute association scores
  names(c(head(x, 5))) 
}))


# Initialize an empty data frame to store changing regulators
changing_regulators_df <- data.frame(Regulators = character(0), F = numeric(0), M = numeric(0))

# Loop through each regulator (row) in viper.markers
for (regulator_name in rownames(viper.markers)) {
  
  # Get the association scores for the current regulator across all clusters
  x <- viper.markers[regulator_name, ]
  
  # Check if the regulator changes from positive to negative or vice versa across all clusters
  if (sum(x > 0) > 0 && sum(x < 0) > 0) {
    # Append the regulator name to the 'Regulators' column
    regulators <- regulator_name
    
    # Append the association scores to the respective columns
    F <- x[1]
    M <- x[2]
    
    # Combine all the data into the data frame
    changing_regulators_df <- rbind(changing_regulators_df, data.frame(
      Regulators = regulators,
      F = F,
      M = M
    ))
  }
}

# Print the resulting data frame
print(changing_regulators_df)
saveRDS(changing_regulators_df$Regulators, "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/Results/changing_regulators_msex_Nvcotfs2_conetwork.rds")


#PA_2plot <- as.matrix(seqad_s[["VIPER"]]@counts[regulators,])
PA_2plot <- as.matrix(seqad_s[["VIPER"]]@counts[changing_regulators_df$Regulators,])
saveRDS(PA_2plot, "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/Results/PA_2plot_New_ntfs11_path.rds")
saveRDS(changing_regulators_df$Regulators, "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/Results/changing_regulators_msex_Nvcotfs2.rds")


#grouping for plots
grouping_column <- paste(seqad_s$Msex, seqad_s$Pathology.group, sep = "-")
grouping_column <- paste(seqad_s$viper.clust, seqad_s$Pathology.group, sep = "-")
grouping_column <- paste(seqad_s$Msex, seqad_s$viper.clust, sep = "-")

# Step 4: Convert the grouping_column to a factor (if it's not already)
grouping_factor <- as.factor(grouping_column)
saveRDS(grouping_factor, "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/UCI/UCI_filtered/Results/grouping_factor_Nvcotfs2.rds")

PA_2plot <- as.matrix(seqad_s[["VIPER"]]@counts[regulators,])
Clusters<-as.factor(as.numeric(seqad_s$viper.clust))
colors<-c("red", "blue") 

my_color_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62")
colors = my_color_palette

# Create PAHeatmap
heatmap_plot <- PAHeatmap(regulators, PA_2plot, Clusters, colors)

# Save PAHeatmap as PDF
ggsave("/lustre/project/hdeng2/nwadiugwu/data/SynapseData/Anderson/Results/ba46_heatmap_plot.pdf", heatmap_plot)


saveRDS(Clusters, "/lustre/project/hdeng2/nwadiugwu/data/SynapseData/Anderson/Results/ba46_Clusters_dianosis.rds")
saveRDS(viper.markers, "viper.markers_pathology.rds")
