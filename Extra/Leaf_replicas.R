#Processing all leaf replicas

leaf.counts_r1<-Read10X("E:/R_images/Leaf_LopezAnido_DevCell/Data/GSM5097889/C1-HCJNCBBXY/filtered_feature_bc_matrix", gene.column = 1)
dim(leaf.counts_r1)
leaf.counts_r2<-Read10X("E:/R_images/Leaf_LopezAnido_DevCell/Data/GSM5097890/C3-HCJNCBBXY/filtered_feature_bc_matrix", gene.column = 1)
dim(leaf.counts_r2)
leaf.counts_r3<-Read10X("E:/R_images/Leaf_LopezAnido_DevCell/Data/GSM5097891/filtered_feature_bc_matrix", gene.column = 1)
dim(leaf.counts_r3)

leaf.dataset_r1 <- CreateSeuratObject(counts = leaf.counts_r1, project = "leaf_r1")

head(leaf.dataset)
dim(leaf.dataset)