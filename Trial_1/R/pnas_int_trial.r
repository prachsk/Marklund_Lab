library(Seurat)
library(future)
library(dplyr)
library(patchwork)
library(ggplot2)
library(remotes)

colon <- readRDS(file = "/Volumes/GoogleDrive/My Drive/Lab/Computations/Trial_1/R/colon.rds")
ileum <- readRDS(file = "/Volumes/GoogleDrive/My Drive/Lab/Computations/Trial_1/R/ileum.rds")
duodenum <- readRDS(file = "/Volumes/GoogleDrive/My Drive/Lab/Computations/Trial_1/R/duodenum.rds")

#all data integration
all.list <- list(colon, ileum, duodenum)
all_features <- SelectIntegrationFeatures(object.list = all.list)

#Perform integration
all.anchors <- FindIntegrationAnchors(object.list = all.list, anchor.features = all_features)

#This command creates an 'integrated' data assay
all.combined <- IntegrateData(anchorset = all.anchors)

# Specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(all.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
all.combined <- ScaleData(all.combined, verbose = FALSE)
all.combined <- RunPCA(all.combined, npcs = 50, verbose = FALSE)
all.combined <- RunUMAP(all.combined, reduction = "pca", dims = 1:30)
all.combined <- FindNeighbors(all.combined, reduction = "pca", dims = 1:30)
all.combined <- FindClusters(all.combined, resolution = 0.5)

p1 <- DimPlot(all.combined, reduction = "umap")
p1

saveRDS(all.combined, file = "/Volumes/GoogleDrive/My Drive/Lab/Computations/Trial_1/R/pnas_int_trial.rds")