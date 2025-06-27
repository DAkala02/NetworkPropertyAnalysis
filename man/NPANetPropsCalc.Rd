\name{NPANetPropsCalc}
\alias{NPANetPropsCalc}
\title{
Parallel Calculation of Spatial Network Properties
}
\description{
Calculates network properties for each gene within each cell using a spatial graph constructed from transcript coordinates and a defined distance threshold. Parallel processing is used to accelerate computation across cell-gene pairs.
}
\usage{
NPANetPropsCalc(
  GroupedFile,
  DistanceThreshold,
  MasterNetworkFile
)
}
\arguments{
  \item{GroupedFile}{
Path to a CSV file containing grouped transcript-level coordinates, with columns {cell_id}, {feature_name}, {x_location}, and {y_location}
  }
  \item{DistanceThreshold}{
Numeric value representing the optimal distance threshold for edge formation in spatial graphs (typically output from {\link{NPADistThresholdDetermination}}).
  }
  \item{MasterNetworkFile}{
File path where the final CSV file of computed network properties will be saved
  }
}
\details{
This function represents the final step in the spatial network pipeline.
It uses parallel computing to process each cell-gene pair in the grouped dataset, constructing a spatial graph by connecting transcripts that lie within the specified {DistanceThreshold}.
For each valid graph, it computes key network statistics such as Mean Degree, Mean Betweenness, Average Path Length, Modularity, Density, and Mean Clustering Coefficient.
Only graphs with at least one edge are analyzed. Results are saved to a CSV file for downstream analysis or visualization.
}
\value{
No return value; function writes a CSV file with the following columns:

\item{CellID}{Identifier of the cell.}
\item{Gene}{Feature name (gene symbol).}
\item{MeanDegree}{Average node degree in the graph.}
\item{MeanBetweenness}{Average betweenness centrality (normalized).}
\item{AvgPathLength}{Mean shortest path length between nodes.}
\item{GraphModularity}{Community modularity based on fast greedy clustering.}
\item{GraphDensity}{Density of the graph (edge-to-node ratio).}
\item{MeanClusteringCoefficient}{Average local clustering coefficient.}
}
\references{
\url{https://igraph.org/r/} \\
\url{https://www.10xgenomics.com/products/xenium}
}
\author{
Daivic B. Akala \email{daivicakala@gmail.com}
}
\note{
This function assumes that the grouped transcript-level input file was created using the earlier step (\code{\link{NPADataAcqOrg}}), and that the optimal distance threshold was selected using \code{\link{NPADistThresholdDetermination}}.
}

\seealso{
\code{\link{NPADataAcqOrg}}, \code{\link{NPADistThresholdDetermination}}, \code{\link[igraph]{graph_from_adjacency_matrix}}
}
\examples{
GroupedFile <- "path/to/GroupedSample.csv"
DistanceThreshold <- #Output from NPADistThresholdDetermination
MasterNetworkFile <- "path/to/NetworkMetrics.csv"

NPANetPropsCalc(
  GroupedFile = GroupedFile,
  DistanceThreshold = DistanceThreshold,
  MasterNetworkFile = MasterNetworkFile
)
}
\keyword{Network}
\keyword{Parallelization}
\keyword{Spatial Transcriptomics}
\concept{Spatial Network Analysis}
\concept{Graph Theory in Biology}
