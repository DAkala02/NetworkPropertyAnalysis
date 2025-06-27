Introduction: NetworkPropertyAnalysis Package

NetworkPropertyAnalysis is an R package designed to process spatial transcriptomic data—particularly from 10x Genomics Xenium—and extract meaningful network-based features at the gene–cell level. 
It includes automated tools to determine an optimal spatial distance threshold, build per-cell/gene graphs, and extract metrics such as degree, betweenness, modularity, and clustering coefficients.
The outputs support downstream analysis including PCA, UMAP, and clustering.

Installation
To install the development version of `NetworkPropertyAnalysis` from GitHub, follow the provided steps:
- Install "devtools"  if not already installed
install.packages("devtools")
- Run the following line in the R Studio console: devtools::install_github("DAkala02/NetworkPropertyAnalysis")
