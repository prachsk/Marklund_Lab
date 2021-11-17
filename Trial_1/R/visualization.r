library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(remotes)

ens <- readRDS(file = "./pnas_loop.rds")
ens_top3 <- head(VariableFeatures(ens), 3)
ens.markers <- FindAllMarkers(ens)
ens.top <- ens.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)

FeaturePlot(object = ens, features = c('Calb2', 'Nos1', 'Chat'),
            pt.size = 0.1,
            sort.cell = T,
            #min.cutoff = 'q1',
            #max.cutoff = 'q90',
            cols = inlmisc::GetColors(99, scheme = 'iridescent', stops = c(0.1,1), bias = 1.0, alpha = 1)
)

VlnPlot(ens, features = c('Elavl4', 'Prph'), pt.size = 0.01)
