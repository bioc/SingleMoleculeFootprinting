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
#' @param verbose Defaults to TRUE
#' 
#' @importFrom parallelDist parDist
#' @importFrom cluster pam silhouette
#' 
#' @export
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
  MethSM_pooled = Reduce(rbind.fill.Matrix, MethSM)
  
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

#' 
#' #' One wrapper function to go from sampleSheet + RegionOfInterest to allelic imbalance according to pam clustering
#' #' 
#' #' @param sampleSheet
#' #' @param RegionOfInterest will be resized in place
#' #' @param TFBSs will be subset for TFBSs overlapping RegionOfInterest (before resize)
#' #' @param params list of parameters to pass to lower functions
#' #' @param deduplicate TRUE or FALSE (default). Useful with samples with high PCR duplication
#' #' @param CytosinesToMask
#' #' @param SMF.plot logical
#' #' @param SMF.plot.params
#' #' 
#' #' @import SingleMoleculeFootprinting
#' #' @import magrittr
#' compute.pam.clustering.Wrapper = function(
#'     sampleSheet, Samples, genome, RegionOfInterest, TFBSs, params,
#'     deduplicate = FALSE, CytosinesToMask = NULL, pool.samples = TRUE, SMF.plot = FALSE){
#'   
#'   # Unpack parameters ####
#'   CallContextMethylation.coverage = params$CallContextMethylation[["coverage"]]
#'   CallContextMethylation.ConvRate.thr = params$CallContextMethylation[["ConvRate.thr"]]
#'   PlottingSNPs = params$PlotSingleSiteSMF[["PlottingSNPs"]]
#'   output.params = params$output.params
#'   ########################
#'   
#'   RegionOfInterest_ext = resize(RegionOfInterest, width = 500, fix = "center")
#'   if(!is.null(TFBSs)){
#'     TFBSs.overlapping.locus = plyranges::filter_by_overlaps(TFBSs[seqnames(TFBSs) == seqnames(RegionOfInterest)], RegionOfInterest_ext)
#'   } else {
#'     TFBSs.overlapping.locus = NULL
#'   }
#'   Experiment = SingleMoleculeFootprinting::DetectExperimentType(Samples = Samples)
#'   
#'   CallContextMethylation(sampleSheet = sampleSheet, sample = Samples,
#'                          genome = genome, RegionOfInterest = RegionOfInterest_ext,
#'                          coverage = CallContextMethylation.coverage, ConvRate.thr = CallContextMethylation.ConvRate.thr,
#'                          returnSM = TRUE, clObj = NULL) -> Methylation
#'   
#'   if (!is.null(CytosinesToMask)){
#'     if(pool.samples){
#'       Methylation <- MaskSNPs(Methylation = Methylation, CytosinesToMask = CytosinesToMask,
#'                               MaskSMmat = TRUE, Experiment = Experiment)
#'     } else {
#'       Methylation <- MaskSNPs2(Methylation = Methylation, CytosinesToMaks = CytosinesToMask,
#'                                MaskSMmat = TRUE, Experiment = Experiment)
#'     }
#'   }
#'   
#'   # NOTE: for now let's just discard endogenous CpG info, as I don't see myself using it in this prj
#'   if (Experiment == "NO"){
#'     Methylation[[1]] = Methylation[[1]]$DGCHN
#'     Methylation[[2]] = lapply(Methylation[[2]], function(x){x$DGCHN})
#'   }
#'   
#'   if(deduplicate){
#'     message("removing PCR duplicates and incomplete molecules")
#'     # source("/g/krebs/barzaghi/Rscripts/R_package/SingleMoleculeFootprinting/R/utilities.r")
#'     lapply(Methylation[[2]], function(x){
#'       incomplete.molecules.idx = as.logical(rowSums(x == 0) > 0)
#'       x = x[!incomplete.molecules.idx,]
#'       duplicated.molecules.idx = as.logical(duplicated(as.matrix(x)))
#'       x = x[!duplicated.molecules.idx,]
#'       return(x)
#'     }) -> Methylation[[2]]
#'     .MethGR = Methylation[[1]]
#'     Methylation[[1]] = SingleMoleculeFootprinting::MethSM.to.MethGR(MethSM =  Methylation[[2]], chromosome = as.character(unique(seqnames(RegionOfInterest_ext))))
#'     overlaps = GenomicRanges::findOverlaps(.MethGR, Methylation[[1]], ignore.strand = TRUE)
#'     strand(Methylation[[1]][subjectHits(overlaps)]) = strand(.MethGR[queryHits(overlaps)])
#'     Methylation[[1]]$GenomicContext = NA
#'     Methylation[[1]][subjectHits(overlaps)]$GenomicContext = .MethGR[queryHits(overlaps)]$GenomicContext
#'     colnames(elementMetadata(Methylation[[1]])) = gsub("_T$", "_Coverage", colnames(elementMetadata(Methylation[[1]])))
#'     colnames(elementMetadata(Methylation[[1]])) = gsub("_M$", "_MethRate", colnames(elementMetadata(Methylation[[1]])))
#'     elementMetadata(Methylation[[1]])[grep("_MethRate$", colnames(elementMetadata(Methylation[[1]])))] =
#'       as.matrix(elementMetadata(Methylation[[1]])[grep("_MethRate$", colnames(elementMetadata(Methylation[[1]])))]) /
#'       as.matrix(elementMetadata(Methylation[[1]])[grep("_Coverage$", colnames(elementMetadata(Methylation[[1]])))])
#'   }
#'   
#'   if (length(Methylation[[1]]) == 0){
#'     message("No cytosines found covering the region")
#'     return()
#'   }
#'   
#'   if (!is.null(PlottingSNPs)){
#'     PlottingSNPs = PlottingSNPs[queryHits(findOverlaps(PlottingSNPs, GRanges(unique(seqnames(Methylation[[1]])), IRanges(min(start(Methylation[[1]])), max(end(Methylation[[1]]))))))]
#'     params$PlotSingleSiteSMF[["PlottingSNPs"]] = PlottingSNPs
#'   }
#'   
#'   if (SMF.plot){
#'     PlotAvgSMF(MethGR = Methylation[[1]], RegionOfInterest = RegionOfInterest_ext, ShowContext = TRUE,
#'                TFBSs = TFBSs.overlapping.locus, SNPs = PlottingSNPs) -> SingleSite.plot
#'   } else {SingleSite.plot = NULL}
#'   
#'   read.origin = lapply(Methylation[[2]], rownames)
#'   MethSM = Reduce(SingleMoleculeFootprinting::rbind.fill.matrix.sparse, Methylation[[2]])
#'   compute.pam.clustering(MethSM = MethSM, TFBSs.overlapping.locus = TFBSs.overlapping.locus, RegionOfInterest = RegionOfInterest, RegionOfInterest_ext = RegionOfInterest_ext,
#'                          params = params, read.origin = read.origin,
#'                          output.params = output.params) -> final.results
#' 
#'   Reduce(rbind,
#'   lapply(seq_along(read.origin), function(i){
#' 
#'     current.sample = names(read.origin)[i]
#'     current.reads = read.origin[[i]]
#'     all.partition.nrs = seq(params$pam[["k"]])
#'     clustering.reads = final.results$pam.res[names(final.results$pam.res) %in% current.reads]
#'     
#'     if(length(clustering.reads) > 0){
#'       
#'       clustering.reads %>%
#'         table() %>%
#'         as.data.frame() %>%
#'         dplyr::rename("partition.nr" = ".", "read.count" = "Freq") %>%
#'         mutate(Sample = current.sample, TFBS.cluster = names(RegionOfInterest)) -> df
#' 
#'       # expand unused levels
#'       if (length(levels(df$partition.nr)) != length(all.partition.nrs)){
#'         missing.partition.nr = all.partition.nrs[!all.partition.nrs %in% df$partition.nr]
#'         df = rbind(
#'           df, 
#'           data.frame(partition.nr = as.factor(missing.partition.nr), read.count = 0, Sample = current.sample, TFBS.cluster = names(RegionOfInterest))
#'         )
#'       }
#'       
#'     } else {
#'       
#'       data.frame(matrix(nrow = 0, ncol = 4, dimnames = list(NULL, c("partition.nr", "read.count", "Sample", "TFBS.cluster")))) -> df
#'     
#'     }
#' 
#'     return(df)
#' 
#'   })) -> cluster.coverage.df
#' 
#'   final.results$cluster.coverage.df = cluster.coverage.df
#'   final.results$read.origin = read.origin
#'   final.results$SingleSite.plot = SingleSite.plot
#'   final.results$RegionOfInterest = RegionOfInterest
#'   
#'   return(final.results)
#'   
#' }
