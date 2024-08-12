#' Convenience wrapper to sort single molecule according to TFBS clusters at multiple sites in the genome
#' 
#' The function starts from a list of single TFBSs, arranges them into clusters, calls methylation at the interested sites and outputs sorted reads
#' 
#' @param sampleSheet QuasR pointer file
#' @param sample samples to use, from the SampleName field of the sampleSheet
#' @param genome BSgenome
#' @param coverage coverage threshold as integer for least number of reads to cover a cytosine for it to be carried over in the analysis. Defaults to 20.
#' @param ConvRate.thr Convesion rate threshold. Double between 0 and 1, defaults to 0.8. To skip this filtering step, set to NULL. For more information, check out the details section.
# #' @param clObj cluster object for parallel processing of multiple samples/RegionsOfInterest. For now only used by qMeth call for bulk methylation. Should be the output of a parallel::makeCluster() call
#' @param CytosinesToMask THIS IS DUCK-TAPE
#' @param TFBSs GRanges object of transcription factor binding sites coordinates
#' @param max_interTF_distance maximum distance between two consecutive TFBSs for them to be grouped in the same window
#' @param max_window_width upper limit to window width. This value should be adjusted according to the user's system as it determines the amount of memory used in the later context methylation call
#' @param min_cluster_width lower limit to window width. Corresponds to the scenario when a window contains a single TFBS.
#' @param bins list of 3 relative bin coordinates. Defaults to list(c(-35,-25), c(-15,15), c(25,35)).
#'             bins[[1]] represents the upstream bin, with coordinates relative to the start of the most upstream TFBS.
#'             bins[[2]] represents all the TFBS bins, with coordinates relative to the center of each TFBS.
#'             bins[[3]] represents the downstream bin, with coordinates relative to the end of the most downstream TFBS.
#' @param sorting_coverage integer. Minimum number of reads covering all sorting bins for sorting to be performed. Defaults to 30.
#' @param cores number of cores to use for parallel processing of multiple Methylation Calling Windows (i.e. groupings of adjecent TFBS clusters)
#' 
#' @importFrom parallel mclapply
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits
#' @importFrom GenomeInfoDb seqlengths
#' 
#' @return list where [[1]] is the TFBSs GRanges object describing coordinates TFBSs used to sort single molecules
#'                    [[2]] is a list of SortedReads nested per TFBS_cluster and sample
#'                    [[3]] is a tibble reporting the count (and frequency) of reads per state, sample and TFBS cluster
#' 
#' @export
#' 
SortReadsBySingleTF_MultiSiteWrapper = function(sampleSheet, sample, genome, coverage = 20, ConvRate.thr = 0.8, # clObj=NULL, ---> parameters passed to CallContextMethylation
                                                 CytosinesToMask = NULL,
                                                 TFBSs,
                                                 max_interTF_distance = 100000, max_window_width = 5000000, min_cluster_width = 600, # ---> parameters passed to Create_MethylationCallingWindows
                                                 sorting_coverage = 30, bins = list(c(-35,-25), c(-15,15), c(25,35)), # ---> parameters passed to SortReadsByTFCluster
                                                 cores = 1
){
  
  names(TFBSs) = paste0("TFBS_", seq(TFBSs))
  
  message("(1) DESIGNING COMMON METHYLATION CALLING WINDOWS FOR ADJACENT CLUSTERS")
  MethylationCallingWindows = Create_MethylationCallingWindows(TFBS_cluster_coordinates = TFBSs,
                                                               max_intercluster_distance = max_interTF_distance,
                                                               max_window_width = max_window_width,
                                                               min_cluster_width = min_cluster_width,
                                                               genomic.seqlenghts = GenomeInfoDb::seqlengths(genome))
  message(paste0(length(MethylationCallingWindows), " METHYLATION CALLING WINDOWS DESIGNED"))
  
  message("(2) CALLING METHYLATION AND SORTING")
  parallel::mclapply(seq_along(MethylationCallingWindows), function(i){
    
    print(i)
    
    CurrentWindow = MethylationCallingWindows[i]
    ExperimentType = suppressMessages(SingleMoleculeFootprinting::DetectExperimentType(Samples = sample))
    
    CallContextMethylation(sampleSheet = sampleSheet,
                           sample = sample,
                           genome = genome,
                           RegionOfInterest = CurrentWindow,
                           coverage = coverage,
                           ConvRate.thr = ConvRate.thr,
                           returnSM = TRUE) -> Methylation
    
    if(length(Methylation[[1]]) == 0){
      return()
    }
    
    if (!is.null(CytosinesToMask)){
      
      message("Masking Cytosines")
      source("/g/krebs/barzaghi/Rscripts/CrappyUtils.R")
      MaskSNPs2(Methylation = Methylation, 
                CytosinesToMaks = CytosinesToMask, 
                MaskSMmat = TRUE, 
                Experiment = ExperimentType) -> Methylation
      
    }
    
    Overlaps = findOverlaps(TFBSs, CurrentWindow)
    TFBSs_to_sort = TFBSs[queryHits(Overlaps)]
    
    lapply(seq_along(TFBSs_to_sort), function(j){
      if(ExperimentType == "NO"){
        Methylation[[2]] = lapply(Methylation[[2]], function(x){x$DGCHN})
      }
      
      SortReadsBySingleTF(MethSM = Methylation[[2]], 
                          TFBS = TFBSs_to_sort[j], 
                          bins = bins, 
                          coverage = sorting_coverage)
    }) -> SortedReads_window
    names(SortedReads_window) = names(TFBSs_to_sort)
    SortedReads_window
    
  }, mc.cores = cores, mc.preschedule = FALSE) -> SortedReads
  SortedReads = unlist(SortedReads, recursive = FALSE)
  
  message("(3) CALCULATE STATE FREQUENCIES")
  Reduce(rbind,
         parallel::mclapply(seq_along(SortedReads), function(i){
           
           StateQuantification_tbl = StateQuantification(SortedReads = SortedReads[[i]], states = NULL)
           StateQuantification_tbl$TFBS_cluster = names(SortedReads[i])
           StateQuantification_tbl
           
         }, mc.cores = cores, mc.preschedule = FALSE)) -> StateFrequency_tbl
  
  return(list(TFBSs, SortedReads, StateFrequency_tbl))
  
}

#' Convenience wrapper to sort single molecule according to TFBS clusters at multiple sites in the genome
#' 
#' The function starts from a list of single TFBSs, arranges them into clusters, calls methylation at the interested sites and outputs sorted reads
#' 
#' @param sampleSheet QuasR pointer file
#' @param sample samples to use, from the SampleName field of the sampleSheet
#' @param genome BSgenome
#' @param coverage coverage threshold as integer for least number of reads to cover a cytosine for it to be carried over in the analysis. Defaults to 20.
#' @param ConvRate.thr Convesion rate threshold. Double between 0 and 1, defaults to 0.8. To skip this filtering step, set to NULL. For more information, check out the details section.
# #' @param clObj cluster object for parallel processing of multiple samples/RegionsOfInterest. For now only used by qMeth call for bulk methylation. Should be the output of a parallel::makeCluster() call
#' @param CytosinesToMask THIS IS DUCK-TAPE
#' @param TFBSs GRanges object of transcription factor binding sites coordinates
#' @param max_intersite_distance maximum allowed distance in base pairs between two TFBS centers for them to be considered part of the same cluster. Defaults to 75.
#' @param min_intersite_distance minimum allowed distance in base pairs between two TFBS centers for them not to be discarded as overlapping. 
#'                               This parameter should be set according to the width of the bins used for later sorting. Defaults to 15.
#' @param max_cluster_size maximum number of TFBSs to be contained in any given cluster. Defaults to 10
#' @param max_intercluster_distance maximum distance between two consecutive TFBS clusters for them to be grouped in the same window
#' @param max_window_width upper limit to window width. This value should be adjusted according to the user's system as it determines the amount of memory used in the later context methylation call
#' @param min_cluster_width lower limit to window width. Corresponds to the scenario when a window contains a single TFBS cluster.
#' @param bins list of 3 relative bin coordinates. Defaults to list(c(-35,-25), c(-7,7), c(25,35)).
#'             bins[[1]] represents the upstream bin, with coordinates relative to the start of the most upstream TFBS.
#'             bins[[2]] represents all the TFBS bins, with coordinates relative to the center of each TFBS.
#'             bins[[3]] represents the downstream bin, with coordinates relative to the end of the most downstream TFBS.
#' @param sorting_coverage integer. Minimum number of reads covering all sorting bins for sorting to be performed. Defaults to 30.
#' @param cores number of cores to use for parallel processing of multiple Methylation Calling Windows (i.e. groupings of adjecent TFBS clusters)
#' 
#' @importFrom parallel mclapply
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits
#' @importFrom GenomeInfoDb seqlengths
#' 
#' @return list where [[1]] is the TFBS_Clusters object describing coordinates and composition of the TFBS clusters used to sort single molecules
#'                    [[2]] is a list of SortedReads nested per TFBS_cluster and sample
#'                    [[3]] is a tibble reporting the count (and frequency) of reads per state, sample and TFBS cluster
#' 
#' @export
#' 
SortReadsByTFCluster_MultiSiteWrapper = function(sampleSheet, sample, genome, coverage = 20, ConvRate.thr = 0.8, # clObj=NULL, ---> parameters passed to CallContextMethylation
                                                 CytosinesToMask = NULL,
                                                 TFBSs, max_intersite_distance = 75, min_intersite_distance = 15, max_cluster_size = 10, max_cluster_width = 300, add.single.TFs = TRUE,  # ---> parameters passed to Arrange_TFBSs_clusters
                                                 max_intercluster_distance = 1e5, max_window_width = 5e6, min_cluster_width = 600, # ---> parameters passed to Create_MethylationCallingWindows
                                                 sorting_coverage = 30, bins = list(c(-35,-25), c(-7,7), c(25,35)), # ---> parameters passed to SortReadsByTFCluster
                                                 cores = 1
                                                 ){
  
  message("(1) CONSTRUCTING TFBS CLUSTERS")
  TFBS_Clusters = Arrange_TFBSs_clusters(TFBSs, 
                                         max_intersite_distance = max_intersite_distance,
                                         min_intersite_distance = min_intersite_distance, 
                                         max_cluster_size = max_cluster_size, 
                                         max_cluster_width = max_cluster_width, 
                                         add.single.TFs = add.single.TFs)
  message(paste0(length(TFBS_Clusters$ClusterCoordinates), " CLUSTERS FOUND"))
  
  message("(2) DESIGNING COMMON METHYLATION CALLING WINDOWS FOR ADJACENT CLUSTERS")
  MethylationCallingWindows = Create_MethylationCallingWindows(TFBS_cluster_coordinates = TFBS_Clusters$ClusterCoordinates,
                                                               max_intercluster_distance = max_intercluster_distance,
                                                               max_window_width = max_window_width,
                                                               min_cluster_width = min_cluster_width,
                                                               genomic.seqlenghts = GenomeInfoDb::seqlengths(genome))
  message(paste0(length(MethylationCallingWindows), " METHYLATION CALLING WINDOWS DESIGNED"))
  
  message("(3) CALLING METHYLATION AND SORTING")
  parallel::mclapply(seq_along(MethylationCallingWindows), function(i){
    
    print(i)
    
    CurrentWindow = MethylationCallingWindows[i]
    ExperimentType = suppressMessages(SingleMoleculeFootprinting::DetectExperimentType(Samples = sample))
    
    CallContextMethylation(sampleSheet = sampleSheet,
                           sample = sample,
                           genome = genome,
                           RegionOfInterest = CurrentWindow,
                           coverage = coverage,
                           ConvRate.thr = ConvRate.thr,
                           returnSM = TRUE) -> Methylation
    
    if(length(Methylation[[1]]) == 0){
      return()
    }
    
    if (!is.null(CytosinesToMask)){
      
      message("Masking Cytosines")
      source("/g/krebs/barzaghi/Rscripts/CrappyUtils.R")
      MaskSNPs2(Methylation = Methylation, 
               CytosinesToMaks = CytosinesToMask, 
               MaskSMmat = TRUE, 
               Experiment = ExperimentType) -> Methylation
      
    }
    
    Overlaps = findOverlaps(TFBS_Clusters$ClusterCoordinates, CurrentWindow)
    Clusters_to_sort = TFBS_Clusters$ClusterComposition[queryHits(Overlaps)]
    
    lapply(seq_along(Clusters_to_sort), function(j){
      if(ExperimentType == "NO"){
        Methylation[[2]] = lapply(Methylation[[2]], function(x){x$DGCHN})
      }
      SortReadsByTFCluster(MethSM = Methylation[[2]],
                           TFBS_cluster = Clusters_to_sort[[j]],
                           bins = bins, 
                           coverage = sorting_coverage)
    }) -> SortedReads_window
    names(SortedReads_window) = names(Clusters_to_sort)
    SortedReads_window
    
  }, mc.cores = cores, mc.preschedule = FALSE) -> SortedReads
  SortedReads = unlist(SortedReads, recursive = FALSE)
  
  message("(4) CALCULATE STATE FREQUENCIES")
  Reduce(rbind,
         parallel::mclapply(seq_along(SortedReads), function(i){
           
           StateQuantification_tbl = StateQuantification(SortedReads = SortedReads[[i]], states = NULL)
           StateQuantification_tbl$TFBS_cluster = names(SortedReads[i])
           StateQuantification_tbl
           
         }, mc.cores = cores, mc.preschedule = FALSE)) -> StateFrequency_tbl
  
  return(list(TFBS_Clusters, SortedReads, StateFrequency_tbl))
  
}

