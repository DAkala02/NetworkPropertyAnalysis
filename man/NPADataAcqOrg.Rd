```r
\name{NPADataAcqOrg}
\alias{NPADataAcqOrg}
\title{
Organize Spatial Transcriptomic Data for Network Property Analysis
}
\description{
Loads raw spatial transcript data from Xenium, organizes it by cell and gene, and saves intermediate files for downstream network analysis
}
\usage{
NPADataAcqOrg(
  TranscriptFile,
  OutputDir,
  GroupedPath,
  GroupedSamplePath,
  SampleSize = 1000,
  Seed = 123
)
}
\arguments{
  \item{TranscriptFile}{
Path to the CSV file containing raw transcript data exported from Xenium. Must include columns: {cell_id}, {transcript_id}, {feature_name}, {x_location}, and {y_location}
  }
  \item{OutputDir}{
Directory where organized intermediate files will be saved; each cell-gene combination will be saved in its own subdirectory
  }
  \item{GroupedPath}{
Path to save the full grouped data (ordered by {cell_id} and {feature_name})
  }
  \item{GroupedSamplePath}{
Path to save a sampled subset of the grouped data used for determining distance thresholds
  }
  \item{SampleSize}{
Number of unique cells to randomly sample for distance threshold analysis; default value is {1000}
  }
  \item{Seed}{
Integer value used to set the random seed for reproducibility; default value is {123}.
  }
}
\details{
This function prepares spatial transcriptomic data for graph-based analysis by grouping transcript
coordinates by {cell_id} and {feature_name}, and organizing them into a nested folder structure.
It also generates a representative sample of cells for threshold calibration in later steps.
Intermediate files are created to facilitate scalable and parallelized downstream computation.
}
\value{
No return value; function is called for its side effects—creating organized CSV files stored in the specified output directory

The following are saved to disk:
  \item{GroupedPath}{CSV containing the full grouped and ordered transcript data.}
  \item{GroupedSamplePath}{CSV containing a subset of sampled cells.}
  \item{OutputDir}{Nested folders containing individual cell-gene CSV files for downstream processing.}
}
\references{
10x Genomics Xenium Platform: \url{https://www.10xgenomics.com/products/xenium}
}
\author{
Daivic B. Akala \email{daivicakala@gmail.com}
}
\note{
This function is designed to be the first step in a three-stage pipeline for computing network properties of spatial transcriptomic data.
}

\seealso{
\code{\link{NPADistThresholdDetermination}}, \code{\link{NPANetPropsCalc}}
}
\examples{
## Example usage:
TranscriptFile <- "path/to/transcripts.csv"
OutputDir <- "path/to/OutputDir"
GroupedPath <- "path/to/GroupedData.csv"
GroupedSamplePath <- "path/to/GroupedSample.csv"

NPADataAcqOrg(
  TranscriptFile = TranscriptFile,
  OutputDir = OutputDir,
  GroupedPath = GroupedPath,
  GroupedSamplePath = GroupedSamplePath,
  SampleSize = 1000,
  Seed = 123
)
}
\keyword{Data}
\keyword{Spatial}
\keyword{Transcriptomics}
\concept{Network Analysis}
```
