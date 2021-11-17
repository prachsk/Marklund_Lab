library(Seurat)
library(dplyr)
library(patchwork)
library(future)
getwd()

folders <- list.files("/Users/pax/Google\ Drive/My\ Drive/Lab/Computations/Trial_1/Data/10X")
ens.list <- list()
result <- "./trial2.rds"
count <- 1

for (files in folders) {
  seurat_data <- Read10X(data.dir = paste0("/Users/pax/Google\ Drive/My\ Drive/Lab/Computations/Trial_1/Data/10X", files))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, min.cells = 3, min.features = 200, project = files)
  seurat_obj <- RenameCells(object = seurat_obj, add.cell.id = paste0("_",count))
  ens.list <- append(ens.list, assign(files, seurat_obj))
  count <- count + 1
}

ens.list <- lapply(X = ens.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = ens.list)

options(future.globals.maxSize = 21000 *1024^2)
plan("multiprocess", workers = 6)
ens.anchors <- FindIntegrationAnchors(object.list = ens.list, anchor.features = features)
ens.combined <- IntegrateData(anchorset = ens.anchors)

DefaultAssay(ens.combined) <- "integrated"

ens.combined <- ScaleData(ens.combined, verbose = FALSE)
ens.combined <- RunPCA(ens.combined, npcs = 50, verbose = FALSE)
ens.combined <- RunUMAP(ens.combined, reduction = "pca", dims = 1:30)
ens.combined <- FindNeighbors(ens.combined, reduction = "pca", dims = 1:30)
ens.combined <- FindClusters(ens.combined, resolution = 0.5)

saveRDS(ens.combined, file = result)