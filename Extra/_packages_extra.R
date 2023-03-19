#Install extra_packages

#install the dependencies for Monocle3
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))
devtools::install_github('cole-trapnell-lab/monocle3')

#install the dependencies for Seurat Wrapper
my_packages_2<- c("R.utils", "spatstat.explore", "spatstat.random")
install.packages(my_packages_2)

devtools::install_github('satijalab/seurat-data')
devtools::install_github("satijalab/seurat-wrappers")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dittoSeq")
