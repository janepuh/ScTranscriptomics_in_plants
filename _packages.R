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

install.packages(my_packages)

devtools::install_github('cole-trapnell-lab/monocle3')
devtools::install_github('renderthis')
install.packages('pdftools')

