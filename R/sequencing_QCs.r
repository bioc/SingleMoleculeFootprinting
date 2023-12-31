#' Conversion rate
#'
#' calculate sequencing library conversion rate on a chromosome of choice
#' @param sampleSheet QuasR sample sheet
#' @param genome BS genome
#' @param chr chromosome to calculate conversion rate on (default: 19)
#' @param clObj cluster object to emply for parallel processing created using the parallel::makeCluster function. Defaults to NULL
#'
#' @importFrom QuasR qMeth
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom BSgenome getSeq
#' @importFrom IRanges resize
#' @importFrom Biostrings vcountPattern
#' @importFrom BiocGenerics grep
#'
#' @return Conversion rate
#'
#' @export
#'
#' @examples
#' Qinput = paste0(tempdir(), "/NRF1Pair_Qinput.txt")
#' library(BSgenome.Mmusculus.UCSC.mm10)
#'
#' if(file.exists(Qinput)){
#'     # DO NOT RUN
#'     # clObj = parallel::makeCluster(5)
#'     # ConversionRatePrecision = ConversionRate(sampleSheet = Qinput, genome = BSgenome.Mmusculus.UCSC.mm10, chr = 19, clObj = clObj)
#'     # parallel::stopCluster(clObj)
#' }
#'
ConversionRate = function(sampleSheet, genome, chr=19, clObj=NULL){

  QuasRprj = GetQuasRprj(sampleSheet, genome)

  # check convertion rate:
  # how many of the non CG/GC cytosines are methylated?
  seq_length = seqlengths(genome)
  chr = tileGenome(seq_length[chr], tilewidth=max(seq_length[chr]), cut.last.tile.in.chrom=TRUE)
  methylation_calls_C = qMeth(QuasRprj, query = chr, mode="allC", reportLevel="C", keepZero = TRUE, clObj = clObj, asGRanges = TRUE, collapseBySample = FALSE)

  seqContext = getSeq(genome, resize(methylation_calls_C, 3, fix='center'))
  GCc = vcountPattern(DNAString("GCN"), seqContext, fixed=FALSE) # take all context for exclusion
  CGc = vcountPattern(DNAString("NCG"), seqContext, fixed=FALSE) # only valid for calling conversion rates
  #SMF
  #non CG non GC
  out_c_met = methylation_calls_C[GCc==0 & CGc==0,]
  # CG_met = methylation_calls_C[!CGc==0,]
  # GC_met = methylation_calls_C[!GCc==0,]
  #non GC context
  # tot.col = grep('R\\d_T',colnames(values(out_c_met)))
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
#' @param sampleSheet QuasR sample sheet
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
#' Qinput = paste0(tempdir(), "/NRF1Pair_Qinput.txt")
#' library(BSgenome.Mmusculus.UCSC.mm10)
#'
#' if(file.exists(Qinput)){
#'     # DO NOT RUN
#'     # clObj = parallel::makeCluster(5)
#'     # BaitRegions = SingleMoleculeFootprintingData::EnrichmentRegions_mm10.rds()
#'     # BaitCaptureEfficiency = BaitCapture(sampleSheet = Qinput, genome = BSgenome.Mmusculus.UCSC.mm10, baits = BaitRegions, clObj=clObj)
#'     # parallel::stopCluster(clObj)
#' }
#'
BaitCapture = function(sampleSheet, genome, baits, clObj=NULL){

  QuasRprj = GetQuasRprj(sampleSheet, genome)

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




# NEED TO DO: following function now doesn't work with external data, and correlation for example data is NA

#' Intersample correlation
#'
#' pair plot of sample correlations
#' @param samples Avg methylation object. Can also be set to "Example" to produce plot using example data of the kind specified by the @param CellType
#' @param context one of "AllCs", "DGCHN", "NWCGW". The first should be chosen for TKO experiments. For experiments carried on WT cells, we recommend checking both the "DGCHN" and "NWCGW" contexts by running this function once per context.
#' @param CellType Cell type to compare your samples to. At the moment, this can be one of "ES", "NP", "TKO".
#'
#' @return Inter-sample correlation plot
#'
#' @export
#'
SampleCorrelation = function(samples, context, CellType){

  # Get methylation data from previous SMF experiments
  # AllC = readRDS(system.file("extdata", "AllCreduced.rds", package = "SingleMoleculeFootprinting", mustWork = TRUE))
  # metMat_ref = readRDS(system.file("extdata", "ReducedRefMat.rds", package = "SingleMoleculeFootprinting", mustWork = TRUE))
  if (!requireNamespace("SingleMoleculeFootprintingData", quietly = TRUE)){
    stop("Install 'SingleMoleculeFootprintingData' to use this function")
    }

  AllC = SingleMoleculeFootprintingData::AllCs.rds()
  metMat_ref = SingleMoleculeFootprintingData::ReferenceMethylation.rds()


  # Pick reference data
  if (length(grep(CellType, colnames(metMat_ref))) > 0){
    refMeth = rowMeans(metMat_ref[,c(grep(CellType, colnames(metMat_ref))[1:2])])
  } else {
    print("Unsupported cell type, skipping comparison to reference dataset")
  }

  # Pick example data
  if (samples == "Example"){
    contextMet = metMat_ref[,c(grep(CellType, colnames(metMat_ref))[3])]
  } else {
    # prepare samples to be a matrix with the right coordinates
    # for (i in seq_along(ncol(contextMet))){
    #   ## Overlap the metCall object with the AllC to define its "absolute" position
    #   ov <- as.matrix(findOverlaps(contextMet[,i], AllC, type = "equal", ignore.strand =F))
    #   ##Use the overlap to fill in the correct positions in the matrix
    #   comparison_mat[,i + 1] [ov[,2]] <- contextMet$met_rate[ov[,1]]
    # }
  }

  # now we don't need this anymore
  remove(metMat_ref)

  # Prepare overall comparison matrix and add samples
  comparison_mat <- matrix(NA, ncol = ncol(contextMet)+1, nrow = length(refMeth))
  comparison_mat[,1] <- refMeth
  comparison_mat[comparison_mat=="NaN"] <- NA
  comparison_mat[,seq_along(ncol(contextMet))] = contextMet
  colnames(comparison_mat) = c(CellType, colnames(contextMet))
  remove(contextMet)

  # Keep track of the rows with only NAs
  NAindex_comparison <- rowSums(is.na(comparison_mat))!=ncol(comparison_mat)

  # Take care of the context
  if (context == "AllCs") {
    CytosinesIndex = rep(TRUE, length(AllC))
  } else if (context == "DGCHN" | context == "NWCGW"){
    CytosinesIndex = AllC$type == context
  } else {
    stop("Unsupported context")
  }

  pairs(comparison_mat[CytosinesIndex & NAindex_comparison,], upper.panel = panel.cor, diag.panel = panel.hist, lower.panel = panel.jet, labels = colnames(comparison_mat))

}
