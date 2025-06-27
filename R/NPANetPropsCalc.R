#' Network Property Analysis for Spatial Transcriptomic Data – Parallel Network Property Calculation
#'
#' Using the determined optimal Distance Threshold, this function computes graph-based network properties
#' for each cell-gene pair and saves the resulting information to a CSV file for additional processing.
#'
#' @param GroupedFile A CSV file of filtered and organized transcript locations with columns:
#' cell_id, feature_name, x_location, y_location.
#' @param DistanceThreshold Numeric distance threshold used to define graph edges (output of previous function).
#' @param MasterNetworkFile File path where the final network metrics CSV will be saved.
#'
#' @return No return value; writes CSV to disk
#'
#' @import data.table
#' @importFrom igraph graph_from_adjacency_matrix degree betweenness mean_distance cluster_fast_greedy modularity edge_density transitivity ecount
#' @importFrom Matrix Matrix
#' @importFrom Rfast Dist
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @export
NPANetPropsCalc <- function(GroupedFile,
                            DistanceThreshold,
                            MasterNetworkFile) {

  # -------------------- Step 1: Load Required Libraries and Initialize Parallel Backend --------------------
  cl <- parallel::makeCluster(parallel::detectCores())
  doParallel::registerDoParallel(cl)

  GroupedData <- data.table::fread(GroupedFile)
  AllCells <- unique(GroupedData$cell_id)
  data.table::setkey(GroupedData, cell_id)

  # -------------------- Step 2: Compute Network Metrics --------------------
  `%dopar%` <- foreach::`%dopar%`
  results <- foreach::foreach(i = seq_along(AllCells),
                              .packages = c("data.table", "igraph", "Matrix", "Rfast")) %dopar% {

                                Cell <- AllCells[i]
                                CellData <- GroupedData[data.table::J(Cell)]
                                if (nrow(CellData) == 0) return(NULL)

                                result_list <- list()

                                for (Gene in unique(CellData$feature_name)) {
                                  GeneData <- CellData[feature_name == Gene]
                                  if (nrow(GeneData) < 1 || nrow(GeneData) > 1000) next

                                  coords <- as.matrix(GeneData[, .(x_location, y_location)])
                                  DistMat <- Rfast::Dist(coords, method = "euclidean", square = TRUE)

                                  AdjMat <- Matrix::Matrix(0, nrow = nrow(GeneData), ncol = nrow(GeneData), sparse = TRUE)
                                  AdjMat[which(DistMat < DistanceThreshold^2 & DistMat > 0, arr.ind = TRUE)] <- 1

                                  Graph <- igraph::graph_from_adjacency_matrix(AdjMat, mode = "undirected")

                                  if (igraph::ecount(Graph) > 0) {
                                    NetStats <- data.table::data.table(
                                      Gene = Gene,
                                      CellID = Cell,
                                      MeanDegree = mean(igraph::degree(Graph)),
                                      MeanBetweenness = mean(igraph::betweenness(Graph, normalized = TRUE)),
                                      AvgPathLength = igraph::mean_distance(Graph, directed = FALSE),
                                      GraphModularity = igraph::modularity(igraph::cluster_fast_greedy(Graph)),
                                      GraphDensity = igraph::edge_density(Graph),
                                      MeanClusteringCoefficient = igraph::transitivity(Graph, type = "average")
                                    )
                                    result_list[[length(result_list) + 1]] <- NetStats
                                  }
                                }

                                data.table::rbindlist(result_list)
                              }

  parallel::stopCluster(cl)

  # -------------------- Step 3: Save Results --------------------
  FinalResults <- data.table::rbindlist(results, use.names = TRUE, fill = TRUE)
  data.table::fwrite(FinalResults, MasterNetworkFile)

  cat("All Network Properties saved into:", MasterNetworkFile, "\n")
}
