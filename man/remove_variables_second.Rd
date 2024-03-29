% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{remove_variables_second}
\alias{remove_variables_second}
\title{Remove variables based on the mean expression levels across samples.}
\usage{
remove_variables_second(sce = NULL, min_meannReads = 0)
}
\arguments{
\item{sce}{A SingleCellExperiment object.}

\item{min_meannReads}{An integer. This function removes variables for which
the mean read counts are less than this value.}
}
\value{
A SingleCellExperiment object.
}
\description{
This function removes variable data such that the mean expression levels
  across samples are less than `min_meannReads`.
}
\examples{
data(pbmc_eg)
pbmc <- remove_variables_second(sce = pbmc_eg, min_meannReads = 0.01)

}
