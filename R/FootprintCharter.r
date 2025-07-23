#' Run FootprintCharter
#' 
#' @param MethSM Single molecule matrix list as returned by CallContextMethylation()
#' @param RegionOfInterest GRanges of coordinates to analyse
#' @param RegionOfInterest_ext RegionOfInterest to be resized, defaults to IRanges::resize(RegionOfInterest, 500, "center")
#' @param TFBSs TFBSs annotation. Used to annotate TF footprints downstream of footprint detection.
#' @param coverage minimum number of molecules required. Defaults to 30
#' @param k number of partitions required. Defaults to 16. Will be dynamically reduced according to minimum number of molecules reqiured (n, see below)
#' @param n minimum number of molecules required per partition
#' @param TF.length vector of two integers for footprint length bounds. Defaults to c(5,75). 
#' @param nucleosome.length vector of two integers for footprint length bounds. Defaults to c(120,1000). 
#' @param cytosine.coverage.thr Cytosine coverage threshold for footprint detection. Individual cytosines will be discarded, not whole footprints. Defaults to 5.
#' @param verbose Defaults to TRUE
#' 
#' @importFrom parallelDist parDist
#' @importFrom cluster pam silhouette
#' @importFrom Seqinfo seqnames
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather
#' @importFrom dplyr right_join select arrange
#' @importFrom magrittr %>%
#' 
#' @export
#' 
#' @examples
#' 
#' Methylation = qs::qread(
#'   system.file("extdata", "Methylation_4.qs", package="SingleMoleculeFootprinting"
#'   ))
#' MethSM = Methylation[[2]]
#' RegionOfInterest = GenomicRanges::GRanges("chr6", IRanges::IRanges(88106000, 88106500))
#' RegionOfInterest = IRanges::resize(RegionOfInterest, 80, "center")
#' 
#' FootprintCharter(
#'   MethSM = MethSM,
#'   RegionOfInterest = RegionOfInterest,
#'   coverage = 30,
#'   k = 16,
#'   n = 5,
#'   TF.length = c(5,75),
#'   nucleosome.length = c(120,1000),
#'   cytosine.coverage.thr = 5,
#'   verbose = TRUE
#'   ) -> FC_results
#'   
FootprintCharter = function(
  MethSM, 
  RegionOfInterest, 
  RegionOfInterest_ext = IRanges::resize(RegionOfInterest, 500, "center"), 
  TFBSs = NULL,
  coverage = 30,
  k = 16,
  n = 5,
  TF.length = c(5,75), 
  nucleosome.length = c(120,1000), 
  cytosine.coverage.thr = 5,
  verbose = TRUE
    ){
  
  if(verbose){message("1. Pooling molecules from all samples")}
  read.origin = lapply(MethSM, rownames)
  MethSM_pooled = Reduce(rbind_fill_sparseMatrix, MethSM)
  
  if(verbose){message("2. Computing sliding windows")}
  MethSM_smoothed = RollingMean(MethSM = MethSM_pooled, RegionOfInterest = RegionOfInterest, verbose = verbose)
  
  if(nrow(MethSM_smoothed) < coverage){
    stop(paste0("The site is covered by less than ", coverage, " continuous reads...quitting"))
  }
  
  if(verbose){message("3. computing distance matrix")}
  distance.matrix = parallelDist::parDist(x = MethSM_smoothed, method = "euclidean", threads = 1)
  
  if(verbose){message("4. Partitioning")}
  partitioned.molecules = cluster::pam(x = distance.matrix, k = k, diss = TRUE, cluster.only = TRUE, keep.diss = FALSE, keep.data = FALSE)
  new.k = k
  while(any(table(partitioned.molecules) < n)){
    if(verbose){message(paste0("partitions too slim detected...retrying with k=", new.k-1))}
    new.k = new.k - 1
    if(new.k == 1){stop('Cannot be clustered given the n required')}
    partitioned.molecules = cluster::pam(x = distance.matrix, k = new.k, diss = TRUE, cluster.only = TRUE, keep.diss = FALSE, keep.data = FALSE)
  }
  
  if(verbose){message("5. Footprints detection")}
  DetectFootprints(
    MethSM = MethSM_pooled, 
    partitioned.molecules = partitioned.molecules, 
    TF.length = TF.length, 
    nucleosome.length = nucleosome.length,
    cytosine.coverage.thr = cytosine.coverage.thr
    ) -> footprints.df
  
  if(!is.null(TFBSs)){
    if(verbose){message("6. Footprints annotation")}
    AnnotateFootprints(
      footprints.df = footprints.df, 
      chromosome = as.character(seqnames(RegionOfInterest)), 
      TFBSs = TFBSs
    ) -> footprints.df
  } else {
    if(verbose){message("6. Footprints annotation (skipping)")}
    footprints.df$seqnames = as.character(seqnames(RegionOfInterest))
    footprints.df$TF = NA
    footprints.df$TF.name = NA
  }
  
  if(verbose){message("7. Footprints aggregation")}
  AggregateFootprints(
    footprints.df = footprints.df
    ) -> footprints.df
  
  if(verbose){message("8. Results wrangling")}
  partitioned.molecules = split(names(partitioned.molecules), unname(partitioned.molecules))
  lapply(read.origin, function(x){
    lapply(partitioned.molecules, function(y){
      x[x%in%y]
      })
    }) -> partitioned.molecules
  
  lapply(partitioned.molecules, lengths) %>%
    data.frame() %>%
    rownames_to_column("partition.nr") %>%
    gather(sample, partition.coverage, -partition.nr) %>%
    right_join(., dplyr::select(footprints.df, -partition.coverage), by = "partition.nr", relationship = "many-to-many") %>%
    arrange(seqnames, start, partition.nr) -> footprints.df
  
  return(list(partitioned.molecules = partitioned.molecules, footprints.df = footprints.df))

}
