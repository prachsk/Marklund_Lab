library(Seurat)
library(future)
library(dplyr)
library(patchwork)
library(ggplot2)
library(remotes)

#Load all data from colon section and create Seurat objects & Visualize normal distribution with histogram plot
colon1.data <- Read10X(data.dir = "/Volumes/GoogleDrive/My Drive/Lab/Computations/Trial_1/Data/10X/GSM4635433_10x-Run0052_Female-Colon-Neurons-Unfixed")
colon1 <- CreateSeuratObject(counts = colon1.data, min.cells = 3 ,min.features = 200, project = "FCNU")
hist(colon1@meta.data$nFeature_RNA, breaks = 500)

colon2.data <- Read10X(data.dir = "/Volumes/GoogleDrive/My Drive/Lab/Computations/Trial_1/Data/10X/GSM4635434_10x-Run0072_Female-Colon-Neurons-Fixed")
colon2 <- CreateSeuratObject(counts = colon2.data, min.cells = 3 ,min.features = 200, project = "FCNF")
hist(colon2@meta.data$nFeature_RNA, breaks = 500)

colon3.data <- Read10X(data.dir = "/Volumes/GoogleDrive/My Drive/Lab/Computations/Trial_1/Data/10X/GSM4635435_10x-Run0060_Male-Colon-Neurons-Fixed")
colon3 <- CreateSeuratObject(counts = colon3.data, min.cells = 3 ,min.features = 200, project = "MCNF")
hist(colon3@meta.data$nFeature_RNA, breaks = 500)

#Violin plots of nFeature_RNA & nCount_RNA from each Seurat object
VlnPlot(colon1, features = c("nFeature_RNA", "nCount_RNA"), pt.size = NULL)
plot1 <- FeatureScatter(colon1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

VlnPlot(colon2, features = c("nFeature_RNA", "nCount_RNA"), pt.size = NULL)
plot2 <- FeatureScatter(colon2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

VlnPlot(colon3, features = c("nFeature_RNA", "nCount_RNA"), pt.size = NULL)
plot3 <- FeatureScatter(colon3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#Set the minimum and maximum cut-off for nFeature_RNA & nCount_RNA in all samples
colon1 <- subset(colon1, subset = nFeature_RNA > 1500 & nFeature_RNA < 9000 & nCount_RNA > 12000 & nCount_RNA < 50000)
colon2 <- subset(colon2, subset = nFeature_RNA > 1000 & nFeature_RNA < 8000 & nCount_RNA > 10000 & nCount_RNA < 60000)
colon3 <- subset(colon3, subset = nFeature_RNA > 1200 & nFeature_RNA < 7500 & nCount_RNA > 10000 & nCount_RNA < 60000)

#Normalize the data with LogNormalize method
colon1 <- NormalizeData(colon1, normalization.method = "LogNormalize", scale.factor = 10000)
colon2 <- NormalizeData(colon2, normalization.method = "LogNormalize", scale.factor = 10000)
colon3 <- NormalizeData(colon3, normalization.method = "LogNormalize", scale.factor = 10000)

#Identification of high variable genes
colon1 <- FindVariableFeatures(colon1, selection.method = "vst", nfeatures = 2000)
colon2 <- FindVariableFeatures(colon2, selection.method = "vst", nfeatures = 2000)
colon3 <- FindVariableFeatures(colon3, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
colon1_top10 <- head(VariableFeatures(colon1), 10)
colon2_top10 <- head(VariableFeatures(colon2), 10)
colon3_top10 <- head(VariableFeatures(colon3), 10)

# plot variable features with and without labels
colon1_plot1 <- VariableFeaturePlot(colon1)
colon1_plot2 <- LabelPoints(plot = colon1_plot1, points = colon1_top10, repel = TRUE)
colon1_plot1 + colon1_plot2

colon2_plot1 <- VariableFeaturePlot(colon2)
colon2_plot2 <- LabelPoints(plot = colon2_plot1, points = colon2_top10, repel = TRUE)
colon2_plot1 + colon2_plot2

colon3_plot1 <- VariableFeaturePlot(colon3)
colon3_plot2 <- LabelPoints(plot = colon3_plot1, points = colon3_top10, repel = TRUE)
colon3_plot1 + colon3_plot2

#Scaling the data
colon1_all.genes <- rownames(colon1)
colon1 <- ScaleData(colon1, features = colon1_all.genes)

colon2_all.genes <- rownames(colon2)
colon2 <- ScaleData(colon2, features = colon2_all.genes)

colon3_all.genes <- rownames(colon3)
colon3 <- ScaleData(colon3, features = colon3_all.genes)

#Run PCA
colon1 <- RunPCA(colon1, features = VariableFeatures(object = colon1))
VizDimLoadings(colon1, dims = 1:2, reduction = "pca")
DimPlot(colon1, reduction = "pca")

colon2 <- RunPCA(colon2, features = VariableFeatures(object = colon2))
VizDimLoadings(colon2, dims = 1:2, reduction = "pca")
DimPlot(colon2, reduction = "pca")

colon3 <- RunPCA(colon3, features = VariableFeatures(object = colon3))
VizDimLoadings(colon3, dims = 1:2, reduction = "pca")
DimPlot(colon3, reduction = "pca")

#Cluster the cells in each dataset
colon1 <- FindNeighbors(colon1, dims = 1:30)
colon1 <- FindClusters(colon1, resolution = 0.5)

colon2 <- FindNeighbors(colon2, dims = 1:30)
colon2 <- FindClusters(colon2, resolution = 0.5)

colon3 <- FindNeighbors(colon3, dims = 1:30)
colon3 <- FindClusters(colon3, resolution = 0.5)

#Run non-linear dimensional reduction (UMAP)
colon1 <- RunUMAP(colon1, dims = 1:30, n.neighbors = 30L, min.dist = 0.5)
DimPlot(colon1, reduction = "umap")

colon2 <- RunUMAP(colon2, dims = 1:30, n.neighbors = 30L, min.dist = 0.5)
DimPlot(colon2, reduction = "umap")

colon3 <- RunUMAP(colon3, dims = 1:30, n.neighbors = 30L, min.dist = 0.5)
DimPlot(colon3, reduction = "umap")

#Find markers for every cluster compared to all remaining cells, report only the positive ones
colon1.markers <- FindAllMarkers(colon1)
colon1.top <- colon1.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
  
colon2.markers <- FindAllMarkers(colon2)
colon2.top <- colon2.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)

colon3.markers <- FindAllMarkers(colon3)
colon3.top <- colon3.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
  
#Colon data integration
colon.list <- list(colon1, colon2, colon3)
colon_features <- SelectIntegrationFeatures(object.list = colon.list)

#Perform integration
colon.anchors <- FindIntegrationAnchors(object.list = colon.list, anchor.features = colon_features)

#This command creates an 'integrated' data assay
colon.combined <- IntegrateData(anchorset = colon.anchors)

# Specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(colon.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
colon.combined <- ScaleData(colon.combined, verbose = FALSE)
colon.combined <- RunPCA(colon.combined, npcs = 50, verbose = FALSE)
colon.combined <- RunUMAP(colon.combined, reduction = "pca", dims = 1:30)
colon.combined <- FindNeighbors(colon.combined, reduction = "pca", dims = 1:30)
colon.combined <- FindClusters(colon.combined, resolution = 0.5)

p1 <- DimPlot(colon.combined, reduction = "umap")
p1

saveRDS(colon.combined, file = "/Volumes/GoogleDrive/My Drive/Lab/Computations/Trial_1/R/colon.rds")

