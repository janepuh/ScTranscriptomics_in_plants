#The script for computer practice 2:
#Seurat: From clusters to marker genes

#Step 0. Preparatory. Loading libraries.
#Use _packages.R to install all the necessary packages prior loading the libraries

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(dittoSeq)

#To make your data calculation reproducible between independent runs, set the seed for the random Number generator
set.seed(100)

#Task 1: Adjusting the number of clusters
#lets repeat the dataset preprocessing from the beginning
path <- 'E:/R_images/Leaf_LopezAnido_DevCell/Data/GSM5097888_Leaf-HVYNNBBXX/filtered_gene_bc_matrices/Arabidopsis/'
leaf.dataset<-Read10X(path, gene.column = 1)
leaf.dataset <- CreateSeuratObject(counts = leaf.dataset, project = "leaf")
leaf.dataset[["percent.mt"]] <- PercentageFeatureSet(leaf.dataset, pattern = "^ATM")
leaf.dataset[["percent.ct"]] <- PercentageFeatureSet(leaf.dataset, pattern = "^ATC")
leaf.dataset <- subset(leaf.dataset, subset = percent.mt <= 20 & percent.ct <= 20 & nCount_RNA >=1000)
leaf.dataset <- SCTransform(leaf.dataset)
leaf.dataset <- RunPCA(leaf.dataset,verbose = FALSE)
leaf.dataset <- RunUMAP(leaf.dataset, dims = 1:50, verbose = FALSE)
leaf.dataset <- FindNeighbors(leaf.dataset, dims = 1:50, verbose = FALSE)

#but this time change resolution parameter in FindClusters function
leaf.dataset_res1 <- FindClusters(leaf.dataset, resolution = 1, verbose = FALSE)
leaf.dataset_res0.2 <- FindClusters(leaf.dataset, resolution = 0.2, verbose = FALSE)
leaf.dataset_res1.5 <- FindClusters(leaf.dataset, resolution = 1.5, verbose = FALSE)

DimPlot(leaf.dataset, label = TRUE, pt.size = 1.5, label.size = 10) + NoLegend()
DimPlot(leaf.dataset_res0.2, label = TRUE, pt.size = 1.5, label.size = 10) + NoLegend()
DimPlot(leaf.dataset_res1.5, label = TRUE, pt.size = 1.5, label.size = 10) + NoLegend()  

#Task 2: Finding marker genes
leaf0.2.markers <- FindAllMarkers(leaf.dataset_res0.2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.csv(leaf0.2.markers, file = "Data/leaf_res0.2_markers.csv")
#leaf0.2.markers<- read.csv(file = "Data/leaf_res0.2_markers.csv", header = TRUE)  

#visualizing marker genes
top2<-leaf0.2.markers %>%
        group_by(cluster) %>%
        slice_max(n = 2, order_by = avg_log2FC)

DoHeatmap(leaf.dataset_res0.2, features = top2$gene) + NoLegend()

#Task 1:Build a heatmap for TOP5 marker genes
