# Install packages

my_packages <- c("tidyverse", 
                 "ggplot2", 
                 "Seurat", 
                 "DoubletDecon", 
                 "xaringan", 
                 "xaringanExtra", 
                 "devtools")

install.packages(my_packages)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dittoSeq")

#install the dependencies for Monocle3
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))
#install the dependencies for Seurat Wrapper
my_packages_2<- c("R.utils", "spatstat.explore", "spatstat.random")
install.packages(my_packages2)

devtools::install_github('satijalab/seurat-data')
devtools::install_github("satijalab/seurat-wrappers")
devtools::install_github('cole-trapnell-lab/monocle3')


devtools::install_github('renderthis')
install.packages('pdftools')

