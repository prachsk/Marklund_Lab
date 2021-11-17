library(Seurat)
library(future)
library(dplyr)
library(patchwork)
library(ggplot2)
library(remotes)

#Load all data from ileum section and create Seurat objects & Visualize normal distribution with histogram plot
ileum1.data <- Read10X(data.dir = "/Volumes/GoogleDrive/My Drive/Lab/Computations/Trial_1/Data/10X/GSM4635439_10x-Run0052_Female-Ileum-Fixed")
ileum1 <- CreateSeuratObject(counts = ileum1.data, min.cells = 3 ,min.features = 200, project = "FIF")
hist(ileum1@meta.data$nFeature_RNA, breaks = 500)

ileum2.data <- Read10X(data.dir = "/Volumes/GoogleDrive/My Drive/Lab/Computations/Trial_1/Data/10X/GSM4635440_10x-Run0052_Female-Ileum-Unfixed")
ileum2 <- CreateSeuratObject(counts = ileum2.data, min.cells = 3 ,min.features = 200, project = "FIU")
hist(ileum2@meta.data$nFeature_RNA, breaks = 500)

ileum3.data <- Read10X(data.dir = "/Volumes/GoogleDrive/My Drive/Lab/Computations/Trial_1/Data/10X/GSM4635441_10x-Run0060_Male-Ileum-Neurons-Fixed")
ileum3 <- CreateSeuratObject(counts = ileum3.data, min.cells = 3 ,min.features = 200, project = "MINF")
hist(ileum3@meta.data$nFeature_RNA, breaks = 500)

ileum4.data <- Read10X(data.dir = "/Volumes/GoogleDrive/My Drive/Lab/Computations/Trial_1/Data/10X/GSM4635442_10x-Run0072-Female-Ileum-Neurons-Fixed")
ileum4 <- CreateSeuratObject(counts = ileum4.data, min.cells = 3 ,min.features = 200, project = "FINF")
hist(ileum4@meta.data$nFeature_RNA, breaks = 500)

#Violin plots of nFeature_RNA & nCount_RNA from each Seurat object
VlnPlot(ileum1, features = c("nFeature_RNA", "nCount_RNA"), pt.size = NULL)
plot1 <- FeatureScatter(ileum1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

VlnPlot(ileum2, features = c("nFeature_RNA", "nCount_RNA"), pt.size = NULL)
plot2 <- FeatureScatter(ileum2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

VlnPlot(ileum3, features = c("nFeature_RNA", "nCount_RNA"), pt.size = NULL)
plot3 <- FeatureScatter(ileum3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

VlnPlot(ileum4, features = c("nFeature_RNA", "nCount_RNA"), pt.size = NULL)
plot4 <- FeatureScatter(ileum4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#Set the minimum and maximum cut-off for nFeature_RNA & nCount_RNA in all samples
ileum1 <- subset(ileum1, subset = nFeature_RNA > 1000 & nFeature_RNA < 8000 & nCount_RNA > 10000 & nCount_RNA < 60000)
ileum2 <- subset(ileum2, subset = nFeature_RNA > 1200 & nFeature_RNA < 7500 & nCount_RNA > 5000 & nCount_RNA < 40000)
ileum3 <- subset(ileum3, subset = nFeature_RNA > 1200 & nFeature_RNA < 8500 & nCount_RNA > 8000 & nCount_RNA < 60000)
ileum4 <- subset(ileum4, subset = nFeature_RNA > 1200 & nFeature_RNA < 7500 & nCount_RNA > 5000 & nCount_RNA < 40000)

#Normalize the data with LogNormalize method
ileum1 <- NormalizeData(ileum1, normalization.method = "LogNormalize", scale.factor = 10000)
ileum2 <- NormalizeData(ileum2, normalization.method = "LogNormalize", scale.factor = 10000)
ileum3 <- NormalizeData(ileum3, normalization.method = "LogNormalize", scale.factor = 10000)
ileum4 <- NormalizeData(ileum4, normalization.method = "LogNormalize", scale.factor = 10000)

#Identification of high variable genes
ileum1 <- FindVariableFeatures(ileum1, selection.method = "vst", nfeatures = 2000)
ileum2 <- FindVariableFeatures(ileum2, selection.method = "vst", nfeatures = 2000)
ileum3 <- FindVariableFeatures(ileum3, selection.method = "vst", nfeatures = 2000)
ileum4 <- FindVariableFeatures(ileum4, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
ileum1_top10 <- head(VariableFeatures(ileum1), 10)
ileum2_top10 <- head(VariableFeatures(ileum2), 10)
ileum3_top10 <- head(VariableFeatures(ileum3), 10)
ileum4_top10 <- head(VariableFeatures(ileum4), 10)

# plot variable features with and without labels
ileum1_plot1 <- VariableFeaturePlot(ileum1)
ileum1_plot2 <- LabelPoints(plot = ileum_plot1, points = ileum_top10, repel = TRUE)
ileum1_plot1 + ileum1_plot2

ileum2_plot1 <- VariableFeaturePlot(ileum2)
ileum2_plot2 <- LabelPoints(plot = ileum2_plot1, points = ileum2_top10, repel = TRUE)
ileum2_plot1 + ileum2_plot2

ileum3_plot1 <- VariableFeaturePlot(ileum3)
ileum3_plot2 <- LabelPoints(plot = ileum3_plot1, points = ileum3_top10, repel = TRUE)
ileum3_plot1 + ileum3_plot2

ileum4_plot1 <- VariableFeaturePlot(ileum4)
ileum4_plot2 <- LabelPoints(plot = ileum4_plot1, points = ileum4_top10, repel = TRUE)
ileum4_plot1 + ileum4_plot2

#Scaling the data
ileum1_all.genes <- rownames(ileum1)
ileum1 <- ScaleData(ileum1, features = ileum1_all.genes)

ileum2_all.genes <- rownames(ileum2)
ileum2 <- ScaleData(ileum2, features = ileum2_all.genes)

ileum3_all.genes <- rownames(ileum3)
ileum3 <- ScaleData(ileum3, features = ileum3_all.genes)

ileum4_all.genes <- rownames(ileum4)
ileum4 <- ScaleData(ileum4, features = ileum4_all.genes)

#Run PCA
ileum1 <- RunPCA(ileum1, features = VariableFeatures(object = ileum1))
VizDimLoadings(ileum1, dims = 1:2, reduction = "pca")
DimPlot(ileum1, reduction = "pca")

ileum2 <- RunPCA(ileum2, features = VariableFeatures(object = ileum2))
VizDimLoadings(ileum2, dims = 1:2, reduction = "pca")
DimPlot(ileum2, reduction = "pca")

ileum3 <- RunPCA(ileum3, features = VariableFeatures(object = ileum3))
VizDimLoadings(ileum3, dims = 1:2, reduction = "pca")
DimPlot(ileum3, reduction = "pca")

ileum4 <- RunPCA(ileum4, features = VariableFeatures(object = ileum4))
VizDimLoadings(ileum4, dims = 1:2, reduction = "pca")
DimPlot(ileum4, reduction = "pca")

#Cluster the cells in each dataset
ileum1 <- FindNeighbors(ileum1, dims = 1:30)
ileum1 <- FindClusters(ileum1, resolution = 0.5)

ileum2 <- FindNeighbors(ileum2, dims = 1:30)
ileum2 <- FindClusters(ileum2, resolution = 0.5)

ileum3 <- FindNeighbors(ileum3, dims = 1:30)
ileum3 <- FindClusters(ileum3, resolution = 0.5)

ileum4 <- FindNeighbors(ileum4, dims = 1:30)
ileum4 <- FindClusters(ileum4, resolution = 0.5)

#Run non-linear dimensional reduction (UMAP)
ileum1 <- RunUMAP(ileum1, dims = 1:30, n.neighbors = 30L, min.dist = 0.5)
DimPlot(ileum1, reduction = "umap")

ileum2 <- RunUMAP(ileum2, dims = 1:30, n.neighbors = 30L, min.dist = 0.5)
DimPlot(ileum2, reduction = "umap")

ileum3 <- RunUMAP(ileum3, dims = 1:30, n.neighbors = 30L, min.dist = 0.5)
DimPlot(ileum3, reduction = "umap")

ileum4 <- RunUMAP(ileum4, dims = 1:30, n.neighbors = 30L, min.dist = 0.5)
DimPlot(ileum4, reduction = "umap")

#Find markers for every cluster compared to all remaining cells, report only the positive ones
ileum1.markers <- FindAllMarkers(ileum1)
ileum1.top <- ileum1.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)

ileum2.markers <- FindAllMarkers(ileum2)
ileum2.top <- ileum2.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)

ileum3.markers <- FindAllMarkers(ileum3)
ileum3.top <- ileum3.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)

ileum4.markers <- FindAllMarkers(ileum4)
ileum4.top <- ileum4.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)

#Ileum data integration
ileum.list <- list(ileum1, ileum2, ileum3, ileum4)
ileum_features <- SelectIntegrationFeatures(object.list = ileum.list)

#Perform integration
ileum.anchors <- FindIntegrationAnchors(object.list = ileum.list, anchor.features = ileum_features)

#This command creates an 'integrated' data assay
ileum.combined <- IntegrateData(anchorset = ileum.anchors)

# Specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(ileum.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
ileum.combined <- ScaleData(ileum.combined, verbose = FALSE)
ileum.combined <- RunPCA(ileum.combined, npcs = 50, verbose = FALSE)
ileum.combined <- RunUMAP(ileum.combined, reduction = "pca", dims = 1:30)
ileum.combined <- FindNeighbors(ileum.combined, reduction = "pca", dims = 1:30)
ileum.combined <- FindClusters(ileum.combined, resolution = 0.5)

p1 <- DimPlot(ileum.combined, reduction = "umap")
p1

saveRDS(ileum.combined, file = "/Volumes/GoogleDrive/My Drive/Lab/Computations/Trial_1/R/ileum.rds")












