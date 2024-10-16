#' Conversion rate
#'
#' calculate sequencing library conversion rate on a chromosome of choice
#' @param sampleFile QuasR sample sheet
#' @param genome BS genome
#' @param chr chromosome to calculate conversion rate on (default: 19)
#' @param cores number of cores for parallel processing. Defaults to 1
#'
#' @importFrom QuasR qMeth
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom parallel makeCluster stopCluster
#' @importFrom BSgenome getSeq
#' @importFrom IRanges resize
#' @importFrom Biostrings vcountPattern
#' @importFrom BiocGenerics grep
#'
#' @export
#'
#' @examples
#'
#' # ConversionRate(sampleFile = sampleFile, 
#' # genome = BSgenome.Mmusculus.UCSC.mm10, chr = 19, cores = 1)
#'
ConversionRate = function(sampleFile, genome, chr=19, cores=1){

  QuasRprj = GetQuasRprj(sampleFile, genome)

  seq_length = seqlengths(genome)
  chr = tileGenome(seq_length[chr], tilewidth=max(seq_length[chr]), cut.last.tile.in.chrom=TRUE)
  cl = makeCluster(cores)
  methylation_calls_C = qMeth(QuasRprj, query = chr, mode="allC", reportLevel="C", keepZero = TRUE, clObj = cl, asGRanges = TRUE, collapseBySample = FALSE)
  stopCluster(cl)

  seqContext = getSeq(genome, trim(resize(methylation_calls_C, 3, fix='center')))
  GCc = vcountPattern(DNAString("GCN"), seqContext, fixed=FALSE) # take all context for exclusion
  CGc = vcountPattern(DNAString("NCG"), seqContext, fixed=FALSE) # only valid for calling conversion rates
  #non CG non GC
  out_c_met = methylation_calls_C[GCc==0 & CGc==0,]
  #non GC context
  tot.col = grep('_T$',colnames(values(out_c_met)))
  met.col = grep('_M$',colnames(values(out_c_met)))
  tot.c = as.matrix(values(out_c_met)[,tot.col])
  met.c = as.matrix(values(out_c_met)[,met.col])
  conv_rate = round((1-(colSums(met.c)/colSums(tot.c)))*100,1)

  return(conv_rate)

}

#' Bait capture efficiency
#'
#' check bait capture efficiency. Expected to be ~70% for mouse genome
#' @param sampleFile QuasR sample sheet
#' @param genome BS genome
#' @param baits GRanges obj of bait coordinates. We provide and example through SingleMoleculeFootprintingData::EnrichmentRegions_mm10.rds()
#' @param clObj cluster object to emply for parallel processing created using the parallel::makeCluster function. Defaults to NULL
#'
#' @import BiocGenerics
#' @importFrom QuasR qCount
#' @importFrom GenomeInfoDb seqlengths
#'
#' @return bait capture efficiency
#'
#' @export
#'
#' @examples
#' 
#' sampleFile = paste0(tempdir(), "/NRF1Pair_Qinput.txt")
#' 
#' if(file.exists(sampleFile)){
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' BaitRegions = SingleMoleculeFootprintingData::EnrichmentRegions_mm10.rds()
#' BaitCapture(sampleFile = sampleFile, genome = BSgenome.Mmusculus.UCSC.mm10, baits = BaitRegions)
#' }
#'
BaitCapture = function(sampleFile, genome, baits, clObj=NULL){

  QuasRprj = GetQuasRprj(sampleFile, genome)

  InBaits=QuasR::qCount(QuasRprj, baits, clObj = clObj)

  seq_length = seqlengths(genome)
  tiles = tileGenome(seq_length, tilewidth = max(seq_length), cut.last.tile.in.chrom=TRUE)
  GW = QuasR::qCount(QuasRprj, tiles, clObj=clObj)

  capture_efficiency = c()
  for(n in seq_along(unique(QuasRprj@alignments$SampleName))){
    capture_efficiency = c(capture_efficiency, sum(InBaits[,n+1]) / sum(GW[,n+1]))
  }

  return(capture_efficiency)

}

#' Plot low coverage methylation rate
#'
#' Inner utility for LowCoverageMethRateDistribution
#'
#' @param Plotting_DF data.frame as returned by GRanges_to_DF function.
#'
#' @import ggplot2
#' @importFrom rlang .data
#'
Plot_LowCoverageMethRate = function(Plotting_DF){

  ggplot(Plotting_DF, aes(x=.data$ExpectedMeth,y=.data$ObservedMeth, group=interaction(.data$Sample,.data$Coverage), color=.data$Sample)) +
    geom_line(aes(linetype=.data$Coverage, size=.data$Coverage)) +
    ylab("Observed Methylation Rate") +
    xlab("Expected Methylation Rate") +
    facet_wrap(~.data$GenomicContext) +
    theme_classic() +
    scale_size_manual("type", values = c(2.5, 1), guide = "none")

}

#' Low Coverage Methylation Rate RMSE
#'
#' Calculate Root mean squared error (RMSE) of methylation rate distribution estimates for low coverage samples
#'
#' @param BinnedMethRate data.frame as returned by GRanges_to_DF function.
#'
#' @importFrom dplyr filter
#'
LowCoverageMethRate_RMSE = function(BinnedMethRate){

  AllSamples = unique(BinnedMethRate$Sample)
  RMSE_DF = data.frame(Sample = c(), RMSE = c())
  for (sample in seq_along(AllSamples)){

    x = filter(BinnedMethRate, Sample == AllSamples[sample])$ExpectedMeth
    y = filter(BinnedMethRate, Sample == AllSamples[sample])$ObservedMeth

    RMSE_DF = rbind(RMSE_DF, data.frame(Sample = AllSamples[sample], RMSE = sqrt(mean((y - x) ** 2))))

  }

  return(RMSE_DF)

}

#' Plot Low Coverage Methylation Rate RMSE
#'
#' Produce barplot of RMSE values calculated for methylation rate distribution estimates of low coverage samples
#'
#' @param RMSE_DF data.frame as returned by the LowCoverageMethRate_RMSE function
#'
#' @import ggplot2
#' @importFrom rlang .data
#'
Plot_LowCoverageMethRate_RMSE = function(RMSE_DF){

  RMSE_DF %>%
    ggplot(aes(.data$Sample, .data$RMSE)) +
    geom_bar(stat = 'identity') +
    ylab("Root mean squared error\n(RMSE)") +
    xlab("") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) -> pl

  return(pl)

}

#' Composite Methylation Rate
#'
#' Monitor methylation rate distribution in a low coverage samples as compared to a high coverage "reference" one.
#' It bins cytosines with similar methylation rates (as observed in the HighCoverage sample) into bins. A single
#' methylation rate value is computed for each bin
#'
#' @param LowCoverage Single GRanges object as returned by CallContextMethylation function run with Coverage parameter set to 1. The object can also contain cytosines from multiple contexts
#' @param LowCoverage_samples Samples to use from the LowCoverage object. Either a string or a vector (for multiple samples).
#' @param HighCoverage Single GRanges object as returned by CallContextMethylation function. The object can also contain cytosines from multiple contexts.
#' @param HighCoverage_samples Single sample to use from HighCoverage. String
#' @param bins The number of bins for which to calculate the "binned" methylation rate. Defaults to 50
#' @param returnDF Whether to return the data.frame used for plotting. Defaults to FALSE
#' @param returnPlot Whether to return the plot. Defaults to TRUE
#' @param RMSE Whether to calculate Mean squared error (RMSE) of methylation rate distribution estimates for low coverage samples. Defaults to TRUE
#' @param return_RMSE_DF  Whether to return a data.frame of computed RMSE values. Defaults to FALSE
#' @param return_RMSE_plot Whether to return a barplot of computed values. Defaults to TRUE
#'
#' @import GenomicRanges
#' @importFrom dplyr mutate
#'
#' @export
#'
#' @examples
#' # I don't have enough example data for this
#' # CompositeMethylationCorrelation(LowCoverage = LowCoverage$DGCHN,
#' #                                 LowCoverage_samples = LowCoverage_Samples,
#' #                                 HighCoverage = HighCoverage$DGCHN,
#' #                                 HighCoverage_samples = HighCoverage_samples[1],
#' #                                 returnDF = FALSE,
#' #                                 returnPlot = TRUE,
#' #                                 RMSE = TRUE,
#' #                                 return_RMSE_DF = FALSE,
#' #                                 return_RMSE_plot = TRUE)
#'
CompositeMethylationCorrelation = function(LowCoverage, LowCoverage_samples, HighCoverage, HighCoverage_samples,
                                           bins = 50, returnDF = FALSE, returnPlot = TRUE,
                                           RMSE = TRUE, return_RMSE_DF = FALSE, return_RMSE_plot = TRUE){

  message("Subsetting GRanges for given samples")
  HighCoverage = SubsetGRangesForSamples(HighCoverage, HighCoverage_samples)
  LowCoverage = SubsetGRangesForSamples(LowCoverage, LowCoverage_samples)

  message("Identifying bins on high coverage sample")
  HighCoverage_MethRate = elementMetadata(HighCoverage)[,grep("_MethRate$", colnames(elementMetadata(HighCoverage)))]
  HighCoverage$Bins = cut(HighCoverage_MethRate, breaks = seq(0,1,1/bins), include.lowest = TRUE)

  message("Assigning bins to low coverage samples")
  Overlaps = findOverlaps(LowCoverage, HighCoverage)
  LowCoverage$Bins = factor(NA, levels = levels(HighCoverage$Bins))
  LowCoverage$Bins[queryHits(Overlaps)] = HighCoverage$Bins[subjectHits(Overlaps)]

  message("Computing binned methylation rates")
  GRanges_objects = list(HighCoverage, LowCoverage)
  Coverages = c("High", "Low")
  BinnedMethRate = Reduce(rbind, lapply(seq_along(GRanges_objects), function(nr){
    mutate(GRanges_to_DF(GRanges_objects[[nr]]), Coverage = Coverages[nr])
  }))

  ReturnList = list()

  if (returnDF){
    ReturnList$MethylationDistribution_DF = BinnedMethRate
  }

  if (returnPlot) {
    Plot = Plot_LowCoverageMethRate(BinnedMethRate)
    ReturnList$MethylationDistribution_plot = Plot
  }

  if (RMSE){

    RMSE_DF = LowCoverageMethRate_RMSE(BinnedMethRate)

    if (return_RMSE_DF){
      ReturnList$RMSE_DF = RMSE_DF
    }

    if (return_RMSE_plot){
      RMSE_plot = Plot_LowCoverageMethRate_RMSE(RMSE_DF)
      ReturnList$RMSE_plot = RMSE_plot
    }

  }

  return(ReturnList)

}