#' Network Property Analysis for Spatial Transcriptomic Data – Optimal Distance Threshold Determination
#'
#' This function processes the developed intermediate files from an output directory and computes network metrics at multiple diameters to determine an optimal threshold using Connectedness and Sparsity metrics.
#'
#' @param CellFile File path to the CSV file that contains the raw cellular data from Xenium
#' @param OutputDir Directory created in the function NPADataAcqOrg containing combined cell+gene CSV files.
#' @param MetricsExcelPath File path to save Excel workbook with network properties for each diameter.
#' @param PlotPath File path which saves the Connectedness vs. Sparsity Plot showing the optimal Distance Threshold
#'
#' @return The estimated optimal Distance Threshold
#'
#' @import data.table
#' @importFrom Matrix Matrix
#' @importFrom igraph graph_from_adjacency_matrix degree vcount
#' @importFrom proxy dist
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook loadWorkbook readWorkbook
#' @importFrom dplyr select filter mutate
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth geom_vline labs theme_minimal
#' @importFrom Rfast Dist
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @export
NPADistThresholdDetermination <- function(CellFile, OutputDir, MetricsExcelPath, PlotPath) {
  # Helper Function 1: Calculate Diameter Range
  get_diameter_range <- function(cell_file) {
    cell_data <- data.table::fread(cell_file)
    diameter <- cell_data$cell_area |>
      na.omit() |>
      (\(x) 2 * sqrt(x / pi))()

    avg_diameter <- mean(diameter, na.rm = TRUE)
    range <- round(min(diameter)):round(max(diameter))
    list(range = range, avg = avg_diameter)
  }

  # Helper Function 2: Process One CSV
  process_csv <- function(file_path, diam, cell_id, gene) {
    df <- tryCatch(data.table::fread(file_path), error = function(e) return(NULL))
    if (is.null(df) || nrow(df) < 2 || !all(c("x_location", "y_location") %in% names(df))) return(NULL)

    dist_mat <- proxy::dist(df[, c("x_location", "y_location")])
    adj <- Matrix::Matrix(0, nrow = nrow(df), ncol = nrow(df), sparse = TRUE)
    adj[which(dist_mat < diam & dist_mat > 0, arr.ind = TRUE)] <- 1

    graph <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
    data.frame(
      CellID = cell_id,
      Gene = gene,
      Diameter = diam,
      MeanDegree = mean(igraph::degree(graph)),
      NumNodes = igraph::vcount(graph),
      NumNA = sum(is.na(adj)),
      NumZero = sum(adj == 0) - igraph::vcount(graph)
    )
  }

  # Helper Function 3: Process One Diameter
  process_diameter <- function(diam) {
    message("Processing Diameter: ", diam)
    SampledCellDirs <- list.dirs(OutputDir, recursive = FALSE, full.names = TRUE)

    metrics <- lapply(SampledCellDirs, function(cell_path) {
      cell_id <- basename(cell_path)
      gene_dirs <- list.dirs(cell_path, recursive = FALSE, full.names = TRUE)

      do.call(rbind, lapply(gene_dirs, function(gene_path) {
        gene <- basename(gene_path)
        csv_files <- list.files(gene_path, pattern = "\\.csv$", full.names = TRUE)
        csv_metrics <- lapply(csv_files, function(file) {
          process_csv(file, diam, cell_id, gene)
        })
        do.call(rbind, csv_metrics)
      }))
    })

    metrics_df <- do.call(rbind, metrics)
    openxlsx::addWorksheet(wb, paste0("D", diam))
    openxlsx::writeData(wb, sheet = paste0("D", diam), x = metrics_df)
    metrics_df
  }

  # Step 1: Compute Diameter Range
  diam_info <- get_diameter_range(CellFile)
  DiameterRange <- diam_info$range
  AverageDiameter <- diam_info$avg

  # Step 2: Create Workbook and Compute Metrics
  wb <- openxlsx::createWorkbook()
  metrics_list <- lapply(DiameterRange, process_diameter)
  openxlsx::saveWorkbook(wb, MetricsExcelPath, overwrite = TRUE)

  # Step 3: Estimate Optimal Threshold
  metric_summary <- sapply(DiameterRange, function(diam) {
    sheetname <- paste0("D", diam)
    df <- tryCatch(openxlsx::readWorkbook(MetricsExcelPath, sheet = sheetname), error = function(e) NULL)
    if (is.null(df) || !all(c("MeanDegree", "NumNodes", "NumNA", "NumZero") %in% names(df))) return(c(NA, NA))

    df <- df |>
      dplyr::mutate(
        Connectedness = MeanDegree / pmax(NumNodes - 1, 1),
        Sparsity = (NumNA + NumZero) / pmax(NumNodes^2, 1)
      )
    c(mean(df$Connectedness, na.rm = TRUE), mean(df$Sparsity, na.rm = TRUE))
  })

  ParameterDF <- data.frame(
    Diameter = DiameterRange,
    Connectedness = metric_summary[1, ],
    Sparsity = metric_summary[2, ]
  ) |> tidyr::drop_na()

  PlotDF <- tidyr::pivot_longer(ParameterDF, cols = c("Connectedness", "Sparsity"), names_to = "Metric", values_to = "Value")
  Plot <- ggplot2::ggplot(PlotDF, ggplot2::aes(x = Diameter, y = Value, color = Metric)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_smooth(method = "loess", se = FALSE, formula = y ~ x, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = AverageDiameter, linetype = "dotdash", color = "black", size = 1) +
    ggplot2::labs(title = "Connectedness vs Sparsity Across Diameters", y = "Metric Value", x = "Diameter") +
    ggplot2::theme_minimal()

  ggplot2::ggsave(PlotPath, Plot)

  LoessConnectedness <- loess(Connectedness ~ Diameter, data = ParameterDF)
  LoessSparsity <- loess(Sparsity ~ Diameter, data = ParameterDF)

  diam_seq <- seq(min(ParameterDF$Diameter), max(ParameterDF$Diameter), by = 0.1)
  y1 <- predict(LoessConnectedness, newdata = data.frame(Diameter = diam_seq))
  y2 <- predict(LoessSparsity, newdata = data.frame(Diameter = diam_seq))
  DistanceThreshold <- round(diam_seq[which.min(abs(y1 - y2))])

  cat(sprintf("Optimal Distance Threshold: %d\n", DistanceThreshold))
  return(DistanceThreshold)
}
