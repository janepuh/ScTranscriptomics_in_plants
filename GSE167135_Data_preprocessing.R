#GSE167135 Data preprocessing

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DropletUtils)

set.seed(42)

#Step 1. Getting data into R
leaf.data<-Read10X("E:/R_images/Leaf_LopezAnido_DevCell/Data/GSM5097888_Leaf-HVYNNBBXX/filtered_gene_bc_matrices/Arabidopsis/",gene.column = 1)
dim(leaf.data)
#there are 31110 genes and 5021 cells in the raw dataset

#Step 2. Making Seurat object
leaf.dataset <- CreateSeuratObject(counts = leaf.data, project = "leaf")

leaf.dataset[["percent.mt"]] <- PercentageFeatureSet(leaf.dataset, pattern = "^ATM")
leaf.dataset[["percent.ct"]] <- PercentageFeatureSet(leaf.dataset, pattern = "^ATC")
##The [[ operator can add columns to object metadata. 
##Check the Seurat object meta.data -> we've got additional variable percent.mt and percent.ct in there calculated for all the cells

leaf.dataset <- subset(leaf.dataset, subset = percent.mt <= 25 & percent.ct <= 25 & nCount_RNA >=500)
dim(leaf.dataset)
#Nothing changed, either because the data was of a good quality or because they already excluded bad cells

VlnPlot(object = leaf.dataset, features = c("nFeature_RNA", "nCount_RNA","percent.ct","percent.mt"), ncol = 5)
#Indeed the authors already filtered out the cells with a low number of counts. Ct and Mt percent looks very good.

#Step 3. Data normalization 
leaf.dataset <- SCTransform(leaf.dataset)
#Note that this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures() functions
#Transformed data will be available in the SCT assay, which is set as the default after running sctransform
#During normalization, we can also remove confounding sources of variation, for example, mitochondrial mapping percentage


#Step 4. Data clustering and visualization
leaf.dataset <- RunPCA(leaf.dataset,verbose = FALSE)
leaf.dataset <- RunUMAP(leaf.dataset, dims = 1:50, verbose = FALSE)

leaf.dataset <- FindNeighbors(leaf.dataset, dims = 1:50, verbose = FALSE)
leaf.dataset <- FindClusters(leaf.dataset, resolution = 1,verbose = FALSE)
DimPlot(leaf.dataset, label = TRUE)
DimPlot(leaf.dataset, label = TRUE, pt.size = 1.5, label.size = 10) + NoLegend()

#Step 5. Gene expression vizualization
FeaturePlot(leaf.dataset,features = 'AT1G77990',order = T,label = F,pt.size = 3)+scale_color_gradientn(colors = c('lightgray','yellow','red','darkred'))

FeaturePlot(leaf.dataset,features = c('AT1G12480', 'AT1G11850', 'AT2G05100', 'AT5G38410'),order = T,label = F,pt.size = 3)

marker_genes<- c('AT5G38420', 'AT1G19150', 'AT2G05100', 'AT5G01530', 'AT4G10340', 'AT1G67090', 'AT5G38410', 'AT3G54890', 'AT1G29920', 'AT3G08940', 'AT3G27690', 'AT5G38430', 'AT1G76570', 'AT4G21750', 'AT2G26330', 'AT3G45640', 'AT2G42840', 'AT5G53210', 'AT5G60880', 'AT2G26330', 'AT1G80080', 'AT1G34245', 'AT3G26744', 'AT1G12860', 'AT2G42840', 'AT5G53210', 'AT5G60880', 'AT2G26330', 'AT5G62230', 'AT5G07180', 'AT1G80080', 'AT1G34245', 'AT1G04110', 'AT4G31805', 'AT3G26744', 'AT1G12860', 'AT3G45640', 'AT3G06120', 'AT5G53210', 'AT2G26330', 'AT5G62230', 'AT5G07180', 'AT1G80080', 'AT1G34245', 'AT2G20875', 'AT1G04110', 'AT4G31805', 'AT3G26744', 'AT1G12860', 'AT3G45640', 'AT3G06120', 'AT2G26330', 'AT5G62230', 'AT5G07180', 'AT2G20875', 'AT1G04110', 'AT4G31805', 'AT3G26744', 'AT1G12860', 'AT1G14350', 'AT3G45640', 'AT3G24140', 'AT3G06120', 'AT2G26330', 'AT5G62230', 'AT5G07180', 'AT2G20875', 'AT1G04110', 'AT3G26744', 'AT1G12860', 'AT1G14350', 'AT3G45640', 'AT3G24140', 'AT2G26330', 'AT5G62230', 'AT5G07180', 'AT1G04110', 'AT1G14350', 'AT5G46240', 'AT3G45640', 'AT3G17300', 'AT5G46240', 'AT2G26330', 'AT3G45640', 'AT4G12970', 'AT2G17950', 'AT2G27250', 'AT2G45190', 'AT2G27250', 'AT4G17710', 'AT1G55580', 'AT2G22850', 'AT3G28730', 'AT4G32880', 'AT5G16560')
# marker_genes<-c('AT1G12480', 'AT4G04890', 'AT5G61480', 'AT5G57350', 'AT4G19200', 'AT1G79430', 'AT1G11850', 'AT1G66400', 'AT5G13930', 'AT5G38420', 'AT1G77990', 'AT5G02380', 'AT5G64700', 'AT5G44160', 'AT3G16330', 'AT1G77110', 'AT2G46570', 'AT2G34710', 'AT2G26580', 'AT5G13930', 'AT1G06360', 'AT1G54040', 'AT4G32880')


DoHeatmap(leaf.dataset, features = marker_genes) + NoLegend()

# Task 1. Clusters annotation. A naive approach. 

leaf.dataset <- RenameIdents(leaf.dataset,
                                                           `0`='Epidermis',
                                                           `1`="Mesophyll",
                                                           `2`='Epidermis',
                                                           `3`="Mesophyll",
                                                           `4`='Phloem',
                                                           `5`='Phloem',
                                                           `6`="Mesophyll",
                                                           `7`='Epidermis',
                                                           `8`='Epidermis',
                                                           `9`='Epidermis',
                                                           `10`='Epidermis',
                                                           `11`='Epidermis',
                                                           `12`='Epidermis',
                                                           `13`='Dividing cells',
                                                           `14`='Xylem',
                                                           `15`='Phloem',
                                                           `16`='Xylem',
                                                           `17`='Sieve element',
                                                           `18`='Phloem')
leaf.dataset$CellTypes<-Idents(leaf.dataset)


Computer practice 2:
#Adjusting the number of clusters 
leaf2.dataset <- FindClusters(leaf.dataset, resolution = 0.5, verbose = FALSE)
DimPlot(leaf.dataset, label = TRUE, pt.size = 1.5, label.size = 10) + NoLegend() + DimPlot(leaf2.dataset, label = TRUE, pt.size = 1.5, label.size = 10) + NoLegend()
