#' Network Property Analysis for Spatial Transcriptomic Data – Data Acquisition and Organization
#'
#' This function loads cell and transcript spatial data from Xenium software and saves organized data into intermediate files within an output directory
#'
#' @param TranscriptFile Describes the file path to the CSV file that contains the raw transcript data from Xenium
#' @param OutputDir Directory created to save intermediate CSV files
#' @param GroupedPath Describes the file path to the CSV file that contains the grouped information from the GroupedData variable
#' @param GroupedSamplePath Describes the file path to the CSV file that contains the grouped information from the SampleData variable
#' @param SampleSize Number of cells to sample for Distance Threshold Determination (Default 1000)
#' @param Seed Random seed value for reproducibility (Default 123)
#'
#' @return No value(s) returned; outputs saved files to appropriate location
#'
#' @import data.table
#' @importFrom dplyr select filter mutate
#' @export
NPADataAcqOrg <- function(TranscriptFile,
                          OutputDir,
                          GroupedPath,
                          GroupedSamplePath,
                          SampleSize = 1000,
                          Seed = 123) {

  # Load Data
  TranscriptData <- data.table::fread(TranscriptFile)
  SelectedData <- data.table::as.data.table(TranscriptData)[, c("cell_id", "transcript_id", "feature_name", "x_location", "y_location"), with = FALSE]
  GroupedData <- SelectedData[order(cell_id, feature_name)]
  if (!is.null(GroupedPath)) {
    data.table::fwrite(GroupedData, GroupedPath, row.names = FALSE)
  }
  data.table::setkey(data.table::setDT(GroupedData), cell_id, feature_name)

  # Sample and Organize
  set.seed(Seed)
  SampledCells <- sample(unique(GroupedData$cell_id), size = SampleSize)
  SampledData <- GroupedData[cell_id %in% SampledCells]

  if (!is.null(GroupedSamplePath)) {
    data.table::fwrite(SampledData, GroupedSamplePath, row.names = FALSE)
  }

  dir.create(OutputDir, recursive = TRUE, showWarnings = FALSE)

  UniqueCells <- unique(SampledData$cell_id)
  nCells <- length(UniqueCells)

  for (i in seq_along(UniqueCells)) {
    Cell <- UniqueCells[i]
    CellDataSubset <- SampledData[cell_id == Cell]
    CellDir <- file.path(OutputDir, paste0("Cell_", Cell))
    dir.create(CellDir, showWarnings = FALSE, recursive = TRUE)

    UniqueGenes <- unique(CellDataSubset$feature_name)

    for (gene in UniqueGenes) {
      GeneData <- CellDataSubset[feature_name == gene]
      if (nrow(GeneData) == 0) next
      GeneDir <- file.path(CellDir, paste0("Gene_", gene))
      dir.create(GeneDir, showWarnings = FALSE, recursive = TRUE)
      OutputFile <- file.path(GeneDir, paste0("Gene_", gene, "_Cell_", Cell, "_Data.csv"))
      data.table::fwrite(GeneData, OutputFile)
    }

    #Progress Reporting
    ProgressPercent <- (i / nCells) * 100
    cat(sprintf("Directory Creation Progress: %.1f%% \r", ProgressPercent))
    flush.console()
  }
}
