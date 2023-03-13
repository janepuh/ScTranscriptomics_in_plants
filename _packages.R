# Install packages

my_packages <- c("tidyverse", 
                 "ggplot2", 
                 "Seurat", 
                 "DoubletDecon", 
                 "xaringan", 
                 "xaringanExtra", 
                 "devtools")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dittoSeq")

#install the dependencies for Monocle3
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))

install.packages(my_packages)

devtools::install_github('cole-trapnell-lab/monocle3')
devtools::install_github('renderthis')
devtools::install_github('cole-trapnell-lab/monocle3', ref="develop")
install.packages('pdftools')

