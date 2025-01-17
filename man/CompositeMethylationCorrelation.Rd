% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sequencing_QCs.r
\name{CompositeMethylationCorrelation}
\alias{CompositeMethylationCorrelation}
\title{Composite Methylation Rate}
\usage{
CompositeMethylationCorrelation(
  LowCoverage,
  LowCoverage_samples,
  HighCoverage,
  HighCoverage_samples,
  bins = 50,
  returnDF = FALSE,
  returnPlot = TRUE,
  RMSE = TRUE,
  return_RMSE_DF = FALSE,
  return_RMSE_plot = TRUE
)
}
\arguments{
\item{LowCoverage}{Single GRanges object as returned by CallContextMethylation function run with Coverage parameter set to 1. The object can also contain cytosines from multiple contexts}

\item{LowCoverage_samples}{Samples to use from the LowCoverage object. Either a string or a vector (for multiple samples).}

\item{HighCoverage}{Single GRanges object as returned by CallContextMethylation function. The object can also contain cytosines from multiple contexts.}

\item{HighCoverage_samples}{Single sample to use from HighCoverage. String}

\item{bins}{The number of bins for which to calculate the "binned" methylation rate. Defaults to 50}

\item{returnDF}{Whether to return the data.frame used for plotting. Defaults to FALSE}

\item{returnPlot}{Whether to return the plot. Defaults to TRUE}

\item{RMSE}{Whether to calculate Mean squared error (RMSE) of methylation rate distribution estimates for low coverage samples. Defaults to TRUE}

\item{return_RMSE_DF}{Whether to return a data.frame of computed RMSE values. Defaults to FALSE}

\item{return_RMSE_plot}{Whether to return a barplot of computed values. Defaults to TRUE}
}
\description{
Monitor methylation rate distribution in a low coverage samples as compared to a high coverage "reference" one.
It bins cytosines with similar methylation rates (as observed in the HighCoverage sample) into bins. A single
methylation rate value is computed for each bin
}
\examples{
# I don't have enough example data for this
# CompositeMethylationCorrelation(LowCoverage = LowCoverage$DGCHN,
#                                 LowCoverage_samples = LowCoverage_Samples,
#                                 HighCoverage = HighCoverage$DGCHN,
#                                 HighCoverage_samples = HighCoverage_samples[1],
#                                 returnDF = FALSE,
#                                 returnPlot = TRUE,
#                                 RMSE = TRUE,
#                                 return_RMSE_DF = FALSE,
#                                 return_RMSE_plot = TRUE)

}
