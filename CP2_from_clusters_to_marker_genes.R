#The script for computer practice 2:
#Seurat: From clusters to marker genes

Computer practice 2:
  #Adjusting the number of clusters 
  leaf2.dataset <- FindClusters(leaf.dataset, resolution = 0.5, verbose = FALSE)
  DimPlot(leaf.dataset, label = TRUE, pt.size = 1.5, label.size = 10) + NoLegend() + DimPlot(leaf2.dataset, label = TRUE, pt.size = 1.5, label.size = 10) + NoLegend()
  