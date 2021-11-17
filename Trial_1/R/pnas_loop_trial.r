library(Seurat)
library(future)
library(dplyr)
library(patchwork)
library(ggplot2)

#list data associated with 10X
folders <- list.files("/Volumes/GoogleDrive/My Drive/Lab/Computations/Trial_1/Data/10X")

#create an empty list
ens.list <- list()
result <- "./pnas_loop.rds"
count <- 1

#read all 10X data in the folder list and create Seurat object associated with the data
#subset data to remove empty droplet
for (files in folders) {
  seurat_data <- Read10X(data.dir = paste0("/Volumes/GoogleDrive/My Drive/Lab/Computations/Trial_1/Data/10X/", files))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, min.cells = 3 ,min.features = 200, project = files)
  seurat_obj <- RenameCells(object = seurat_obj, add.cell.id = paste0("_",count))
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 1200 & nFeature_RNA < 7500 & nCount_RNA > 8000 & nCount_RNA < 65000)
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  ens.list <- append(ens.list, assign(files, seurat_obj)) 
  count <- count + 1 
}

options(future.globals.maxSize = 3000 *1024^2)
plan("multicore", workers = 6)

#pnas ens data integration
ens_features <- SelectIntegrationFeatures(object.list = ens.list)

#Perform integration
ens.anchors <- FindIntegrationAnchors(object.list = ens.list, anchor.features = ens_features)

#This command creates an 'integrated' data assay
ens.combined <- IntegrateData(anchorset = ens.anchors)

# Specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(ens.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
ens.combined <- ScaleData(ens.combined, verbose = FALSE)
ens.combined <- RunPCA(ens.combined, npcs = 50, verbose = FALSE)
ens.combined <- RunUMAP(ens.combined, reduction = "pca", dims = 1:30)
ens.combined <- FindNeighbors(ens.combined, reduction = "pca", dims = 1:30)
ens.combined <- FindClusters(ens.combined, resolution = 0.5)

p1 <- DimPlot(ens.combined, reduction = "umap")
p1

saveRDS(ens.combined, file = result)