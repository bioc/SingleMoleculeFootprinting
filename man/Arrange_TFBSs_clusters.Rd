% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methylation_calling_utils.r
\name{Arrange_TFBSs_clusters}
\alias{Arrange_TFBSs_clusters}
\title{Convenience function to arrange a list of given TFBSs into clusters}
\usage{
Arrange_TFBSs_clusters(
  TFBSs,
  max_intersite_distance = 75,
  min_intersite_distance = 15,
  max_cluster_size = 6,
  max_cluster_width = 300,
  add.single.TFs = TRUE
)
}
\arguments{
\item{TFBSs}{GRanges object of TFBSs}

\item{max_intersite_distance}{maximum allowed distance in base pairs between two TFBS centers for them to be considered part of the same cluster. Defaults to 75.}

\item{min_intersite_distance}{minimum allowed distance in base pairs between two TFBS centers for them not to be discarded as overlapping. 
This parameter should be set according to the width of the bins used for later sorting. Defaults to 15.}

\item{max_cluster_size}{maximum number of TFBSs to be contained in any given cluster. Defaults to 6}

\item{max_cluster_width}{maximum width of TFBS clusters in bps. Defaults to 300}

\item{add.single.TFs}{Whether to add the TFs not used to create TFBS.clusters to the list for sorting. Defaults to TRUE}
}
\value{
list with two elements: ClusterCoordinates (GRanges object of clusters coordinates) and ClusterComposition (GRangesList of sites for each cluster)
}
\description{
For each TFBS, the genomic neighborhood defined by max_cluster_width will be scanned for adjacent TFBSs. 
The hits will be filtered for min_intersite_distance where, in case of overlapping TFBSs, the second TFBS will be arbitrarily dropped.
These TFBSs plus the central "anchoring" one will define a TFBS cluster. 
This approach implies that the same TFBS can be employed to design multiple clusters in a sliding-window fashion.
}
\examples{

KLF4s = qs::qread(system.file("extdata", "KLF4_chr19.qs", package="SingleMoleculeFootprinting"))
Arrange_TFBSs_clusters(KLF4s)

}
