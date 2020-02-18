# https://github.com/agitter/single-cell-pseudotime

install.packages("VGAM")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("slingshot")
BiocManager::install("destiny")
install.packages('DeLorean')

BiocManager::install("monocle")
# for monocle 3
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor'))
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')


install.packages("SCORPIUS")
install.packages("SLICER")
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("zji90/TSCAN")

install.packages('lle')
devtools::install_github("kieranrcampbell/pseudogp")
install.packages("princurve", type = "source")

BiocManager::install("TSCAN")

install.packages("energy")
install.packages("foreach")
