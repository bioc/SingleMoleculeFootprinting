#' Collect bulk SMF data for later composite plotting
#'
#' @param sampleFile QuasR sampleFile
#' @param samples vector of unique sample names corresponding to the SampleName field from the sampleFile
#' @param genome BSgenome
#' @param TFBSs GRanges object of TF binding sites to collect info for. We reccommend employing 50 to 200 TFBSs.
#' @param window window size to collect methylation information for
#' @param coverage coverage threshold as integer for least number of reads to cover a cytosine for it to be carried over in the analysis. Defaults to 20.
#' @param ConvRate.thr Convesion rate threshold. Double between 0 and 1, defaults to NULL For more information, check out the details section
#' @param cores number of cores to use
#'
#' @import GenomicRanges
#' @import dplyr
#' @importFrom tidyr separate
#' @importFrom parallel makeCluster stopCluster
#' @importFrom rlang .data
#'
#' @return data.frame of bulk SMF info ready for plotting
#'
#' @export
#' 
#' @examples
#' 
#' sampleFile = NULL
#' if(!is.null(sampleFile)){
#' CollectCompositeData(
#' sampleFile = sampleFile, 
#' samples = samples, 
#' genome = BSgenome.Mmusculus.UCSC.mm10, 
#' TFBSs = TopMotifs, 
#' window = 1000, 
#' coverage = 20, 
#' ConvRate.thr = NULL, 
#' cores = 16
#' ) -> CompositeData
#' }
#'
CollectCompositeData = function(sampleFile, samples, genome, TFBSs, window, coverage=20, ConvRate.thr = NULL, cores=1){

  RegionOfInterest = IRanges::resize(TFBSs, window, "center")
  
  clObj = parallel::makeCluster(cores)
  CallContextMethylation(
    sampleFile = sampleFile, 
    samples = samples, 
    genome = genome, 
    RegionOfInterest = RegionOfInterest, 
    coverage = coverage, 
    ConvRate.thr = ConvRate.thr,
    returnSM = FALSE, 
    clObj = clObj, 
    verbose = FALSE
  ) -> Methylation
  parallel::stopCluster(clObj)
  
  if(DetectExperimentType(samples) == "NO"){Methylation = Methylation$DGCHN}
  if(length(Methylation) == 0){return(NULL)}
  
  findOverlaps(Methylation, RegionOfInterest, ignore.strand = TRUE) -> Overlaps
  if(length(Overlaps) == 0){return(Methylation)}
  
  # This also expands cytosines overlapping with more than 1 TFBS    
  Methylation_expanded = Methylation[queryHits(Overlaps)]
  Methylation_expanded$TFBS_center = NA  
  Methylation_expanded$TFBS_center = start(IRanges::resize(RegionOfInterest[subjectHits(Overlaps)], width = 1, fix = 'center'))
  Methylation_expanded$TFBS_index = subjectHits(Overlaps)
  
  message(paste0(round(sum(table(Methylation_expanded$TFBS_index) > 1)/length(TFBSs)*100), 
                 "% of sites are found covered at least at 1 cytosine"))
  
  Methylation_expanded %>%
    data.frame() %>%
    mutate(TFBS_strandedness = as.character(strand(RegionOfInterest))[.data$TFBS_index]) %>%
    mutate(RelStart = ifelse(.data$TFBS_strandedness == "+", start - .data$TFBS_center, -(start - .data$TFBS_center))) %>%
    dplyr::select(-.data$TFBS_strandedness) %>%
    gather(.data$Measure, .data$Value, grep("_Coverage$|_MethRate$", colnames(elementMetadata(Methylation_expanded)), value = TRUE)) %>%
    mutate(Sample = gsub("_MethRate$|_Coverage$", "", .data$Measure), Measure = gsub("^.*_", "", .data$Measure)) %>%
    spread(.data$Measure, .data$Value) %>%
    mutate(SMF = 1 - .data$MethRate) %>%
    dplyr::select(-.data$MethRate) -> CompositeData
  
  return(CompositeData)

}


#' Plot composite SMF data
#'
#' Will use geom_point with <= 5000 points, geom_hex otherwise 
#' 
#' @param CompositeData the output of the CollectCompositeData function
#' @param span the span parameter to pass to geom_smooth
#' @param TF string of TF name to use for plot title
#' 
#' @import ggplot2
#' @import viridis
#' @importFrom ggpointdensity geom_pointdensity
#' @importFrom rlang .data
#' 
#' @export
#' 
#' @examples
#' 
#' # CompositePlot(CompositeData = CompositeData, span = 0.1, TF = "Rest")
#' 
CompositePlot = function(CompositeData, span=0.1, TF){
  
  if(nrow(CompositeData) <= 5000){
    point_density = FALSE
  } else {
    point_density = TRUE
  }
  
  CompositeData %>%
    ggplot(aes(.data$RelStart, .data$SMF)) +
    {if(!point_density){geom_point(aes(color = .data$Sample), alpha=0.5)}} +
    {if(point_density){geom_pointdensity(adjust = 0.1, size = .5)}} +
    geom_smooth(se = TRUE, method = "loess", span = span) +
    viridis::scale_fill_viridis(option = "inferno") +
    ylim(c(0,1)) +
    xlab("Coord relative to TFBS center") +
    ggtitle(paste0("Top ", length(unique(CompositeData$TFBS_index)), " covered ", TF, " sites")) +
    facet_wrap(~Sample) +
    {if(point_density){viridis::scale_color_viridis(trans = "log10")}} +
    theme_classic() 
  
}







