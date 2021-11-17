library(Seurat)
library(future)
library(dplyr)
library(patchwork)
library(ggplot2)
library(remotes)

#Load all data from duodenum section and create Seurat objects & Visualize normal distribution with histogram plot
duodenum1.data <- Read10X(data.dir = "/Volumes/GoogleDrive/My Drive/Lab/Computations/Trial_1/Data/10X/GSM4635445_10x-Run0060_Male-Doudenum-Neurons_Fixed")
duodenum1 <- CreateSeuratObject(counts = duodenum1.data, min.cells = 3 ,min.features = 200, project = "FCNU")
hist(duodenum1@meta.data$nFeature_RNA, breaks = 500)

duodenum2.data <- Read10X(data.dir = "/Volumes/GoogleDrive/My Drive/Lab/Computations/Trial_1/Data/10X/GSM4635446_10x-Run0072_Female-Duodenum-Neurons-Fixed")
duodenum2 <- CreateSeuratObject(counts = duodenum2.data, min.cells = 3 ,min.features = 200, project = "FCNF")
hist(duodenum2@meta.data$nFeature_RNA, breaks = 500)

#Violin plots of nFeature_RNA & nCount_RNA from each Seurat object
VlnPlot(duodenum1, features = c("nFeature_RNA", "nCount_RNA"), pt.size = NULL)
plot1 <- FeatureScatter(duodenum1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

VlnPlot(duodenum2, features = c("nFeature_RNA", "nCount_RNA"), pt.size = NULL)
plot2 <- FeatureScatter(duodenum2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#Set the minimum and maximum cut-off for nFeature_RNA & nCount_RNA in all samples
duodenum1 <- subset(duodenum1, subset = nFeature_RNA > 1500 & nFeature_RNA < 9000 & nCount_RNA > 10000 & nCount_RNA < 80000)
duodenum2 <- subset(duodenum2, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & nCount_RNA > 8000 & nCount_RNA < 50000)

#Normalize the data with LogNormalize method
duodenum1 <- NormalizeData(duodenum1, normalization.method = "LogNormalize", scale.factor = 10000)
duodenum2 <- NormalizeData(duodenum2, normalization.method = "LogNormalize", scale.factor = 10000)

#Identification of high variable genes
duodenum1 <- FindVariableFeatures(duodenum1, selection.method = "vst", nfeatures = 2000)
duodenum2 <- FindVariableFeatures(duodenum2, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
duodenum1_top10 <- head(VariableFeatures(duodenum1), 10)
duodenum2_top10 <- head(VariableFeatures(duodenum2), 10)

# plot variable features with and without labels
duodenum1_plot1 <- VariableFeaturePlot(duodenum1)
duodenum1_plot2 <- LabelPoints(plot = duodenum1_plot1, points = duodenum1_top10, repel = TRUE)
duodenum1_plot1 + duodenum1_plot2

duodenum2_plot1 <- VariableFeaturePlot(duodenum2)
duodenum2_plot2 <- LabelPoints(plot = duodenum2_plot1, points = duodenum2_top10, repel = TRUE)
duodenum2_plot1 + duodenum2_plot2

#Scaling the data
duodenum1_all.genes <- rownames(duodenum1)
duodenum1 <- ScaleData(duodenum1, features = duodenum1_all.genes)

duodenum2_all.genes <- rownames(duodenum2)
duodenum2 <- ScaleData(duodenum2, features = duodenum2_all.genes)

#Run PCA
duodenum1 <- RunPCA(duodenum1, features = VariableFeatures(object = duodenum1))
VizDimLoadings(duodenum1, dims = 1:2, reduction = "pca")
DimPlot(duodenum1, reduction = "pca")

duodenum2 <- RunPCA(duodenum2, features = VariableFeatures(object = duodenum2))
VizDimLoadings(duodenum2, dims = 1:2, reduction = "pca")
DimPlot(duodenum2, reduction = "pca")

#Cluster the cells in each dataset
duodenum1 <- FindNeighbors(duodenum1, dims = 1:30)
duodenum1 <- FindClusters(duodenum1, resolution = 0.5)

duodenum2 <- FindNeighbors(duodenum2, dims = 1:30)
duodenum2 <- FindClusters(duodenum2, resolution = 0.5)

#Run non-linear dimensional reduction (UMAP)
duodenum1 <- RunUMAP(duodenum1, dims = 1:30, n.neighbors = 30L, min.dist = 0.5)
DimPlot(duodenum1, reduction = "umap")

duodenum2 <- RunUMAP(duodenum2, dims = 1:30, n.neighbors = 30L, min.dist = 0.5)
DimPlot(duodenum2, reduction = "umap")

#Find markers for every cluster compared to all remaining cells, report only the positive ones
duodenum1.markers <- FindAllMarkers(duodenum1)
duodenum1.top <- duodenum1.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)

duodenum2.markers <- FindAllMarkers(duodenum2)
duodenum2.top <- duodenum2.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)

#duodenum data integration
duodenum.list <- list(duodenum1, duodenum2)
duodenum_features <- SelectIntegrationFeatures(object.list = duodenum.list)

#Perform integration
duodenum.anchors <- FindIntegrationAnchors(object.list = duodenum.list, anchor.features = duodenum_features)

#This command creates an 'integrated' data assay
duodenum.combined <- IntegrateData(anchorset = duodenum.anchors)

# Specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(duodenum.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
duodenum.combined <- ScaleData(duodenum.combined, verbose = FALSE)
duodenum.combined <- RunPCA(duodenum.combined, npcs = 50, verbose = FALSE)
duodenum.combined <- RunUMAP(duodenum.combined, reduction = "pca", dims = 1:30)
duodenum.combined <- FindNeighbors(duodenum.combined, reduction = "pca", dims = 1:30)
duodenum.combined <- FindClusters(duodenum.combined, resolution = 0.5)

p1 <- DimPlot(duodenum.combined, reduction = "umap")
p1

saveRDS(duodenum.combined, file = "/Volumes/GoogleDrive/My Drive/Lab/Computations/Trial_1/R/duodenum.rds")