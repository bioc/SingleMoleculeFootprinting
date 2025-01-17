#' Get QuasRprj
#'
#' @param sampleFile QuasR pointer file
#' @param genome BSgenome
#'
#' @import QuasR
#'
#' @export
#'
#' @examples
#'
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' CacheDir <- ExperimentHub::getExperimentHubOption(arg = "CACHE")
#' sampleFile = paste0(CacheDir, "/NRF1Pair_sampleFile.txt")
#' QuasRprj = GetQuasRprj(sampleFile, BSgenome.Mmusculus.UCSC.mm10)
#'
GetQuasRprj = function(sampleFile, genome){

  QuasRprj=QuasR::qAlign(sampleFile=sampleFile,
                        genome=genome@pkgname,
                        projectName = "prj",
                        paired="fr",
                        aligner = "Rbowtie",
                        bisulfite="undir")
  QuasRprj@aligner = "Rbowtie"

  return(QuasRprj)

}

#' Get Single Molecule methylation matrix
#'
#' Used internally as the first step in CallContextMethylation
#'
#' @param QuasRprj QuasR project object as returned by calling the QuasR function qAling on previously aligned data
#' @param range GenimocRange representing the genomic region of interest
#' @param sample One of the sample names as reported in the SampleName field of the QuasR sampleFile provided to qAlign. N.b. all the files with the passed sample name will be used to call methylation
#'
#' @import QuasR
#' @importFrom Matrix sparseMatrix rsparsematrix
#' @importFrom IRanges isEmpty
#'
#' @return List of single molecule methylation matrixes (all Cytosines), one per sample
#'
#' @export
#'
#' @examples
#' 
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' library(IRanges)
#' library(GenomicRanges)
#' 
#' CacheDir <- ExperimentHub::getExperimentHubOption(arg = "CACHE")
#' sampleFile = paste0(CacheDir, "/NRF1Pair_sampleFile.txt")
#' sample = suppressMessages(readr::read_delim(sampleFile, delim = "\t")[[2]])
#' QuasRprj = GetQuasRprj(sampleFile, BSgenome.Mmusculus.UCSC.mm10)
#' range = GRanges("chr6", IRanges(88106000, 88106500))
#' 
#' GetSingleMolMethMat(QuasRprj, range, sample)
#'
GetSingleMolMethMat = function(QuasRprj,range,sample){

  QuasRprj_sample = QuasRprj[unlist(lapply(unique(sample), function(s){grep(s, QuasRprj@alignments$SampleName)}))]
  Cs=QuasR::qMeth(QuasRprj_sample, query=range, mode="allC",reportLevel="alignment", collapseBySample = TRUE)

  MethMatrix_list = lapply(seq_along(Cs), function(i){
    
    if (any(isEmpty(Cs[[i]]))){
      message(paste0("No single molecule methylation info found in the given range for sample ", names(Cs)[i]))
      meth_matrix = Matrix::rsparsematrix(nrow = 0, ncol = 0, density = 0)
      meth_matrix
    } else {
      ReadID = factor(Cs[[i]]$aid)
      Coordinate = factor(Cs[[i]]$Cid)
      meth_matrix = Matrix::sparseMatrix(i = as.numeric(ReadID),
                                         j = as.numeric(Coordinate),
                                         x = Cs[[i]]$meth + 1,
                                         dimnames = list(levels(ReadID), as.character(levels(Coordinate))))
      meth_matrix@factors = list(strand = unlist(lapply(split(Cs[[i]]$strand, Cs[[i]]$Cid), unique), use.names = FALSE))
      meth_matrix
    }
      
  })
  
  names(MethMatrix_list) = sort(sample)

  return(MethMatrix_list)
}

#' Compute MethGR from MethSM
#'
#' @param MethSM internal CallContextMethylation
#' @param chromosome string
#'
MethSM.to.MethGR = function(MethSM, chromosome){
  
  lapply(seq_along(MethSM), function(i){
    
    if(any(dim(MethSM[[i]])) > 0){
      MethGR = GRanges(seqnames = chromosome, IRanges(start = as.integer(colnames(MethSM[[i]])), width = 1), strand = MethSM[[i]]@factors$strand)
      MethGR$`_T` = Matrix::colSums(MethSM[[i]]>0)
      MethGR$`_M` = Matrix::colSums(MethSM[[i]]==2)
      colnames(elementMetadata(MethGR)) = paste0(names(MethSM)[i], colnames(elementMetadata(MethGR)))
    } else {
      MethGR = GRanges("mock", IRanges(1,2))
      MethGR$`_T` = NA
      MethGR$`_M` = NA
      colnames(elementMetadata(MethGR)) = paste0(names(MethSM)[i], colnames(elementMetadata(MethGR)))
      MethGR = MethGR[0]
    }
    
    return(MethGR)
    
  }) -> MethGR.list
  
  MethGR = Reduce(full.join.granges, MethGR.list)
  tmp = as.matrix(elementMetadata(MethGR))
  tmp[is.na(tmp)] = 0 
  elementMetadata(MethGR) = tmp
  MethGR = sort(MethGR)
  
  return(MethGR)
  
}

#' Calculate reads conversion rate
#'
#' @param MethSM as comes out of the func GetSingleMolMethMat
#' @param chr Chromosome, MethSM doesn't carry this info
#' @param genome BSgenome
#' @param thr Double between 0 and 1. Threshold below which to filter reads.
#'
#' @import GenomicRanges
#' @import Biostrings
#' @importFrom IRanges IRanges
#'
#' @return Filtered MethSM
#'
#' @export
#'
#' @examples
#'
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' MethSM = qs::qread(system.file("extdata", "Methylation_3.qs", 
#' package="SingleMoleculeFootprinting"))[[2]]$SMF_MM_TKO_DE_
#' FilterByConversionRate(MethSM, chr = "chr19", 
#' genome = BSgenome.Mmusculus.UCSC.mm10, thr = 0.8)
#'
FilterByConversionRate = function(MethSM, chr, genome, thr){

  CytosineRanges = GRanges(rep(chr, ncol(MethSM)),IRanges(as.numeric(colnames(MethSM)),width = 1))
  GenomicContext = Biostrings::getSeq(genome, resize(CytosineRanges,3,fix='center'))
  IsInContext = Biostrings::vcountPattern('GC',GenomicContext)==1 | vcountPattern('CG',GenomicContext)==1

  # filter based on conversion
  ConvRate = 1 - (rowMeans_drop0(MethSM[,!IsInContext, drop=FALSE]) - 1)
  if (nrow(MethSM) != length(ConvRate)){
    stop("Discrepancy during FilterByConversionRate execution")
  }
  message(paste0(100 - round((sum(ConvRate>thr, na.rm = TRUE)/length(ConvRate))*100, digits = 2), "% of reads found with conversion rate below ", thr))
  if(sum(is.na(ConvRate))>0){
    message(paste0(sum(is.na(ConvRate)), " reads have no cytosines out of context. For these, conversion rate cannot be computed and they will be discarded"))
  }
  ConvRate[is.na(ConvRate)] = 0
  FilteredSM = MethSM[ConvRate>thr,, drop=FALSE]

  return(FilteredSM)

}

#' Detect type of experiment
#'
#' @param Samples SampleNames field from QuasR sampleFile
#'
#' @export
#'
#' @examples
#'
#' CacheDir = ExperimentHub::getExperimentHubOption(arg = "CACHE")
#' sampleFile = paste0(CacheDir, "/NRF1Pair_sampleFile.txt")
#' samples = suppressMessages(unique(readr::read_delim(sampleFile, delim = "\t")[[2]]))
#' DetectExperimentType(samples)
#'
DetectExperimentType = function(Samples){

  if(length(grep("_NO_", Samples)) > 0){
    ExpType = "NO"
  }else if(length(grep("_SS_",Samples)) > 0){
    ExpType = "SS"
  }else if(length(grep("_DE_",Samples)) > 0){
    ExpType = "DE"
  }else{stop("Sample name error. No _McvPI_ , _SsssI_ or _McvPISsssI_")}

  message(paste0("Detected experiment type: ", ExpType))
  return(ExpType)

}

#' Filter Cytosines in context
#'
#' @param MethGR Granges obj of average methylation
#' @param genome BSgenome
#' @param context Context of interest (e.g. "GC", "CG",..)
#'
#' @import GenomicRanges
#' @import Biostrings
#'
#' @return filtered Granges obj
#'
#' @export
#'
#' @examples
#'
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' MethGR = qs::qread(system.file("extdata", "Methylation_3.qs", 
#' package="SingleMoleculeFootprinting"))[[1]]
#'
#' FilterContextCytosines(MethGR, BSgenome.Mmusculus.UCSC.mm10, "NGCNN")
#'
FilterContextCytosines = function(MethGR, genome, context){

  CytosineRanges = GRanges(seqnames(MethGR), ranges(MethGR), strand(MethGR)) # performing operation without metadata to make it lighter
  seqinfo(CytosineRanges) = seqinfo(MethGR)

  NewWidth = ifelse(nchar(context) %% 2 == 1, nchar(context), nchar(context) + 1) # Make context nchar() odd so we can fix to 'center'
  GenomicContext = Biostrings::getSeq(genome, suppressWarnings(trim(IRanges::resize(CytosineRanges, width = NewWidth, fix = "center"))))# I checked: it is strand-aware | , as.character = FALSE I don't know what was this here for

  # Fix the truncated sites
  trimmed <- which(!width(GenomicContext) == NewWidth)
  for(k in trimmed){
    if(as.character(strand(CytosineRanges[k])) == '+'){
      GenomicContext[k] <- paste0(c(rep('N', 5-width(GenomicContext[k])), as.character(GenomicContext[k])), collapse="")
    } else if(as.character(strand(CytosineRanges[k])) == '-'){
      GenomicContext[k] <- paste0(c(as.character(GenomicContext[k]), rep('N', 5-width(GenomicContext[k]))), collapse="")
    }
  }

  FixedPattern = !(length(grep("A|C|G|T", strsplit(gsub("N", "", context), "")[[1]], invert = TRUE)) > 0) # check whether we need to use fixed or not in the next function...fixed==TRUE is a tiny bit faster
  IsInContext = Biostrings::vcountPattern(context, subject = GenomicContext, fixed = FixedPattern) > 0

  MethGR_InContext = MethGR[IsInContext]
  MethGR_InContext$GenomicContext = rep(context, length(MethGR_InContext))

  return(MethGR_InContext)

}

#' Collapse strands
#'
#' @param MethGR Granges obj of average methylation
#' @param context "GC" or "HCG". Broad because indicates just the directionality of collapse.
#'
#' @import GenomicRanges
#'
#' @return MethGR with collapsed strands (everything turned to - strand)
#' 
#' @export
#' 
#' @examples
#' 
#' # CollapseStrands(MethGR, "GC")
#' 
#'
CollapseStrands = function(MethGR, context){

  if(length(MethGR) == 0){
    return(MethGR)
  }
  
  if(context != "GC" & context != "HCG"){
    warning("Unrecognized context, please pass one of GC and HCG")
  }
  
  if (length(unique(MethGR$GenomicContext)) > 1){
    stop("There should be only one value for GenomicContext here...more will cause problems")
  }

  # find the - stranded cytosines and make them +
  MethGR_minus = MethGR[strand(MethGR) == "-"]
  start(MethGR_minus) = start(MethGR_minus) + ifelse(grepl("HCG|CG", context), -1, +1)
  end(MethGR_minus) = start(MethGR_minus)
  strand(MethGR_minus) = "+"

  # Extract the + Cs
  MethGR_plus = MethGR[strand(MethGR) == "+"]

  # Sum the counts
  ov = findOverlaps(MethGR_minus, MethGR_plus)
  values(MethGR_plus[subjectHits(ov)])[,-length(values(MethGR))] = as.matrix(values(MethGR_plus[subjectHits(ov)])[,-length(values(MethGR_plus)),drop=FALSE]) + as.matrix(values(MethGR_minus[queryHits(ov)])[,-length(values(MethGR_minus)),drop=FALSE])
  if(length(queryHits(ov)) == 0){ # N.b. Negative indexing with integer(0) returns and empty object!!
    CollapsedMethGR = sort(c(MethGR_plus, MethGR_minus), by = ~ seqnames + start + end)
  } else {
    CollapsedMethGR = sort(c(MethGR_plus, MethGR_minus[-queryHits(ov)]), by = ~ seqnames + start + end)
  }
  

  return(CollapsedMethGR)

}

#' Collapse strands in single molecule matrix
#'
#' The idea here is that (regardless of context) if a C is on the - strand, calling getSeq on that coord (N.b. unstranded, that's the important bit) will give a "G', a "C" if it's a + strand.
#'
#' @param MethSM Single molecule matrix
#' @param context "GC" or "CG". Broad because indicates just the directionality of collapse.
#' @param genome BSgenome
#' @param chr Chromosome, MethSM doesn't carry this info
#'
#' @import GenomicRanges
#' @import Biostrings
#' @importFrom IRanges IRanges
#' @importFrom Matrix rowSums
#'
#' @return Strand collapsed MethSM
#' 
#' @export
#' 
#' @examples
#' 
#' # CollapseStrandsSM(MethSM, "GC", BSgenome.Mmusculus.UCSC.mm10, "chr19")
#' 
#'
CollapseStrandsSM = function(MethSM, context, genome, chr){

  if(any(MethSM@Dim == 0)){
    return(MethSM)
  }

  CytosineRanges = GRanges(chr,IRanges(as.numeric(colnames(MethSM)),width = 1))
  GenomicContext = Biostrings::getSeq(genome, CytosineRanges)
  IsMinusStrand = GenomicContext=="G" # no need to select for GC/CG context --> MethSM passed already subset

  # Separate MethSM into - and + strands
  MethSM_minus = MethSM[Matrix::rowSums(MethSM[,IsMinusStrand, drop=FALSE]) > 0, IsMinusStrand, drop=FALSE]
  MethSM_plus = MethSM[Matrix::rowSums(MethSM[,!IsMinusStrand, drop=FALSE]) > 0, !IsMinusStrand, drop=FALSE]
  NrMinReads = dim(MethSM_minus)[1]
  message(paste0(ifelse(is.null(NrMinReads), 0, NrMinReads), " reads found mapping to the - strand, collapsing to +"))

  # Turn - into +
  offset = ifelse(grepl("GC", context), +1, -1) # the opposite if I was to turn + into -
  colnames(MethSM_minus) = as.character(as.numeric(colnames(MethSM_minus)) + offset)

  # Merge matrixes
  StrandCollapsedSM = rbind.fill.Matrix(x = MethSM_minus, y = MethSM_plus)

  return(StrandCollapsedSM)

}

#' Filter Cs for coverage
#'
#' @param MethGR Granges obj of average methylation
#' @param thr converage threshold
#'
#' @import GenomicRanges
#'
#' @return filtered MethGR
#' 
#' @export
#'
CoverageFilter = function(MethGR, thr){

  Tcounts = as.matrix(elementMetadata(MethGR)[grep("_T$", colnames(elementMetadata(MethGR)))])
  Mcounts = as.matrix(elementMetadata(MethGR)[grep("_M$", colnames(elementMetadata(MethGR)))])

  # filter for coverage
  CoverageMask = Tcounts < thr
  Tcounts[CoverageMask] = NA
  Mcounts[CoverageMask] = NA

  # Calculate methylation rate
  MethylationMatrix = Mcounts/Tcounts

  # bind the GRanges with the scores
  elementMetadata(MethGR)[grep("_T$", colnames(elementMetadata(MethGR)))] = Tcounts
  elementMetadata(MethGR)[grep("_M$", colnames(elementMetadata(MethGR)))] = MethylationMatrix
  colnames(elementMetadata(MethGR)) = gsub("_T$", "_Coverage", colnames(elementMetadata(MethGR)))
  colnames(elementMetadata(MethGR)) = gsub("_M$", "_MethRate", colnames(elementMetadata(MethGR)))

  # do the actual filtering
  if(length(MethGR) > 0){
    MethGRFiltered = MethGR[!(rowSums(is.na(elementMetadata(MethGR)[,-length(elementMetadata(MethGR))])) == (length(elementMetadata(MethGR))-1))]
    return(MethGRFiltered)
  } else {
    return(MethGR)
  }
    
}

#' Call Context Methylation
#'
#' Can deal with multiple samples
#'
#' @param sampleFile QuasR pointer file
#' @param samples vector of unique sample names corresponding to the SampleName field from the sampleFile
#' @param genome BSgenome
#' @param RegionOfInterest GenimocRange representing the genomic region of interest
#' @param coverage coverage threshold as integer for least number of reads to cover a cytosine for it to be carried over in the analysis. Defaults to 20.
#' @param ConvRate.thr Convesion rate threshold. Double between 0 and 1, defaults to NULL. To skip this filtering step, set to NULL. For more information, check out the details section.
#' @param returnSM whether to return the single molecule matrix, defaults to TRUE
#' @param clObj cluster object for parallel processing of multiple samples. For now only used by qMeth call for bulk methylation. Should be the output of a parallel::makeCluster() call
#' @param verbose whether to print out messages while executing. Defaults to FALSE
#'
#' @import QuasR
#' @import GenomicRanges
#' @import BiocGenerics
#'
#' @return List with two Granges objects: average methylation call (GRanges) and single molecule methylation call (matrix)
#' 
#' @details The ConvRate.thr argument should be used with care as it could create biases (e.g. when only one C out of context is present) while generally only marginally cleaning up the data. 
#'
#' @export
#' 
#' @examples
#'
#' sampleFile = NULL
#' if(!is.null(sampleFile)){
#' Methylation <- CallContextMethylation(
#'   sampleFile = sampleFile, 
#'   samples = samples, 
#'   genome = BSgenome.Mmusculus.UCSC.mm10, 
#'   RegionOfInterest = RegionOfInterest,
#'   coverage = 20, 
#'   returnSM = TRUE,
#'   ConvRate.thr = NULL,
#'   clObj = NULL
#' )
#' }
#'
CallContextMethylation = function(sampleFile, samples, genome, RegionOfInterest, coverage = 20, ConvRate.thr = NULL, returnSM = TRUE, clObj=NULL, verbose = FALSE){

  if(verbose){message("Setting QuasR project")}
  QuasRprj = GetQuasRprj(sampleFile, genome)

  if(verbose){message("Calling methylation at all Cytosines")}
  QuasRprj_sample = QuasRprj[unlist(lapply(unique(samples), function(s){grep(s, QuasRprj@alignments$SampleName)}))]
  if (returnSM){
    if(verbose){message("Extracting single molecule matrix")}
    MethSM = GetSingleMolMethMat(QuasRprj, RegionOfInterest, samples)
    if (!is.null(ConvRate.thr)){
      MethSM_ConvRateFiltered = lapply(seq_along(unique(samples)), function(i){
        FilterByConversionRate(MethSM[[i]], chr = seqnames(RegionOfInterest), genome = genome, thr = ConvRate.thr)})
      MethSM = MethSM_ConvRateFiltered
    }
    MethGR = MethSM.to.MethGR(MethSM = MethSM, chromosome = unique(seqnames(RegionOfInterest)))
  } else {
    MethGR = QuasR::qMeth(QuasRprj_sample, mode="allC", query = RegionOfInterest, collapseBySample = TRUE, keepZero = TRUE, clObj = clObj) %>% sort()
  }
  
  if(verbose){message("checking if RegionOfInterest contains information at all")}
  CoverageCols = grep("_T$", colnames(elementMetadata(MethGR)))
  lapply(CoverageCols, function(i){
    sum(elementMetadata(MethGR)[,i]) > 0
  }) -> Samples_covered
  if (all(!unlist(Samples_covered))){
    if(verbose){message("No bulk methylation info found for the given RegionOfInterest for any of the samples")}
    MethGR = MethGR[elementMetadata(MethGR)[,1]!=0]
    colnames(elementMetadata(MethGR)) = gsub("_T$", "_Coverage", colnames(elementMetadata(MethGR)))
    colnames(elementMetadata(MethGR)) = gsub("_M$", "_MethRate", colnames(elementMetadata(MethGR)))
    if (returnSM){
      MethSM = list(Matrix::rsparsematrix(nrow=0,ncol=0,density = 0))
      return(list(MethGR, MethSM))
    } else {
      return(MethGR)
    }
  }

  if(verbose){message("Discard immediately the cytosines not covered in any sample")}
  CsToKeep = rowSums(as.matrix(elementMetadata(MethGR)[,CoverageCols])) > 0
  MethGR = MethGR[CsToKeep]

  if(verbose){message("Subsetting Cytosines by permissive genomic context (GC, HCG)")}
  # Here we use a permissive context: needed for the strand collapsing
  ContextFilteredMethGR = list(GC = FilterContextCytosines(MethGR, genome, "GC"),
                               CG = FilterContextCytosines(MethGR, genome, "HCG"))
  if (returnSM){
    ContextFilteredMethSM = lapply(seq_along(MethSM),
                                   function(n){lapply(seq_along(ContextFilteredMethGR),
                                                      function(i){MethSM[[n]][,colnames(MethSM[[n]]) %in% as.character(start(ContextFilteredMethGR[[i]])), drop=FALSE]})})
  }
  
  if(verbose){message("Collapsing strands")}
  StrandCollapsedMethGR = list(GC = CollapseStrands(MethGR = ContextFilteredMethGR[[1]], context = "GC"),
                               CG = CollapseStrands(MethGR = ContextFilteredMethGR[[2]], context = "HCG"))
  if (returnSM){
    StrandCollapsedMethSM = lapply(seq_along(ContextFilteredMethSM),
                                   function(n){
                                     list(GC = CollapseStrandsSM(ContextFilteredMethSM[[n]][[1]], context = "GC", genome = genome, chr = as.character(seqnames(RegionOfInterest))),
                                          CG = CollapseStrandsSM(ContextFilteredMethSM[[n]][[2]], context = "HCG", genome = genome, chr = as.character(seqnames(RegionOfInterest))))})
  }
  
  if(verbose){message("Filtering Cs for coverage")}
  CoverageFilteredMethGR = list(GC = CoverageFilter(MethGR = StrandCollapsedMethGR[[1]], thr = coverage),
                                CG = CoverageFilter(MethGR = StrandCollapsedMethGR[[2]], thr = coverage))
  
  if (returnSM){
    CoverageFilteredMethSM = lapply(seq_along(StrandCollapsedMethSM),
                                    function(n){lapply(seq_along(CoverageFilteredMethGR),
                                                       function(i){
                                                         SampleCoverageColumn = grep("_Coverage$", colnames(elementMetadata(CoverageFilteredMethGR[[i]])))[n]
                                                         if(is.na(SampleCoverageColumn)){SampleCoverageColumn = NULL}
                                                         CsCoveredEnough = as.character(start(CoverageFilteredMethGR[[i]]))[
                                                           !is.na(data.frame(elementMetadata(CoverageFilteredMethGR[[i]])[,SampleCoverageColumn]))]
                                                         x = StrandCollapsedMethSM[[n]][[i]][,colnames(StrandCollapsedMethSM[[n]][[i]]) %in% CsCoveredEnough, drop=FALSE]
                                                         if (any(dim(x) == 0)){x = Matrix::rsparsematrix(nrow=0,ncol=0,density = 0)}
                                                         x
                                                       })})
  }
  
  # Determining strict context based on ExpType
  ExpType = DetectExperimentType(samples)
  if (ExpType == "NO"){
    ExpType_contexts = c("DGCHN", "-", "NWCGW")
  } else if (ExpType == "SS"){
    ExpType_contexts = c("-", "-", "CG")
  } else if (ExpType == "DE"){
    ExpType_contexts = c("GCH", "GCG", "HCG")
  }

  if(verbose){message(paste0("Subsetting Cytosines by strict genomic context (", paste(ExpType_contexts, collapse = ", "),") based on the detected experiment type: ", ExpType))}
  ContextFilteredMethGR_strict = list(FilterContextCytosines(CoverageFilteredMethGR[[1]], genome, ExpType_contexts[1]),
                                      FilterContextCytosines(CoverageFilteredMethGR[[1]], genome, ExpType_contexts[2]),
                                      FilterContextCytosines(CoverageFilteredMethGR[[2]], genome, ExpType_contexts[3]))
  names(ContextFilteredMethGR_strict) = ExpType_contexts
  if (returnSM){
    ContextFilteredMethSM_strict = lapply(seq_along(CoverageFilteredMethSM), function(sample){ 
      list(CoverageFilteredMethSM[[sample]][[1]][,colnames(CoverageFilteredMethSM[[sample]][[1]]) %in% as.character(start(ContextFilteredMethGR_strict[[1]])), drop=FALSE],
           CoverageFilteredMethSM[[sample]][[1]][,colnames(CoverageFilteredMethSM[[sample]][[1]]) %in% as.character(start(ContextFilteredMethGR_strict[[2]])), drop=FALSE],
           CoverageFilteredMethSM[[sample]][[2]][,colnames(CoverageFilteredMethSM[[sample]][[2]]) %in% as.character(start(ContextFilteredMethGR_strict[[3]])), drop=FALSE])
    })
  }

  if (ExpType == "DE"){
    if(verbose){message("Merging matrixes")}
    MergedGR = sort(Reduce(append, ContextFilteredMethGR_strict), by = ~ seqnames + start + end)
    NAMES = unique(gsub("_Coverage$", "", grep("_Coverage$", colnames(elementMetadata(MergedGR)), value=TRUE)))
    if (returnSM){
      MergedSM = lapply(seq_along(ContextFilteredMethSM_strict), function(n){
        Reduce(cbind.fill.Matrix, ContextFilteredMethSM_strict[[n]])
        })
      }
  } else {
    if(verbose){message("Returning non merged matrixes")}
    MockContext = names(ContextFilteredMethGR_strict) != "-"
    MergedGR = ContextFilteredMethGR_strict[MockContext]
    NAMES = unique(gsub("_Coverage$", "", grep("_Coverage$", colnames(elementMetadata(MergedGR$DGCHN)), value=TRUE)))
    if (returnSM){
      MergedSM = lapply(ContextFilteredMethSM_strict, function(x){x[MockContext]})
      for (i in seq_along(MergedSM)){
        names(MergedSM[[i]]) = ExpType_contexts[MockContext]
        }
      }
  }

  if (returnSM){
    names(MergedSM) = NAMES
    return(list(MergedGR, MergedSM))
  } else {
    return(MergedGR)
  }

}
