#####
# LowCoverage_MethRate_SampleCorrelation utils
#

#' Subset Granges for given samples
#'
#' Inner utility for LowCoverageMethRateDistribution
#'
#' @param GRanges_obj GRanges object as returned by CallContextMethylation function
#' @param Samples vector of sample names as they appear in the SampleName field of the QuasR sampleFile
#'
#' @import GenomicRanges
#' @importFrom dplyr as_tibble
#'
SubsetGRangesForSamples = function(GRanges_obj, Samples){
  
  Cols2Discard = grep(paste0(c(Samples, "^GenomicContext$"), collapse = "|"), colnames(elementMetadata(GRanges_obj)), invert = TRUE)
  elementMetadata(GRanges_obj)[,Cols2Discard] = NULL
  All_NA_rows = is.na(rowMeans(as_tibble(elementMetadata(GRanges_obj)[,-ncol(elementMetadata(GRanges_obj))]), na.rm = TRUE))
  GRanges_obj_subset = GRanges_obj[!All_NA_rows]
  
  return(GRanges_obj_subset)
  
}

#' Manipulate GRanges into data.frame
#'
#' Inner utility for LowCoverageMethRateDistribution
#'
#' @param GRanges_obj GRanges object as returned by CallContextMethylation function
#'
#' @import dplyr
#' @importFrom tidyr spread gather extract
#' @importFrom stats na.omit
#' @importFrom rlang .data
#'
GRanges_to_DF = function(GRanges_obj){
  
  GRanges_obj %>%
    as_tibble() %>%
    dplyr::select(matches("^seqnames$|^start$|_MethRate$|_Coverage$|^GenomicContext$|^Bins$")) %>%
    gather(Sample, Score, -.data$seqnames, -.data$start, -.data$GenomicContext, -.data$Bins) %>%
    na.omit() %>%
    extract(.data$Sample, c("Sample", "Measurement"), regex = "(.*)_(Coverage$|MethRate|)") %>%
    spread(.data$Measurement, .data$Score) %>%
    mutate(Methylated = .data$Coverage*.data$MethRate) %>%
    group_by(.data$Sample, .data$GenomicContext, .data$Bins) %>%
    summarise(TotCoverage = sum(.data$Coverage), TotMethylated = sum(.data$Methylated), ObservedMeth = .data$TotMethylated/.data$TotCoverage) %>%
    ungroup() %>%
    group_by(.data$Sample, .data$GenomicContext) %>%
    mutate(ExpectedMeth = seq(1,100,100/dplyr::n())/100) %>%
    ungroup() %>%
    select(-.data$TotCoverage, -.data$TotMethylated, -.data$Bins) -> DF
  
  return(DF)
  
}

#####
# HighCoverage_MethRate_SampleCorrelation utils
#

#' Utility for HighCoverage_MethRate_SampleCorrelation
#' 
#' @param ... data for lower pairs panel
#'
#' @import graphics
#' 
#' @export
#'
panel.jet <- function(...) {
  jet.colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  smoothScatter(..., nrpoints=0, add=TRUE, colramp=jet.colors) }

#' Utility for HighCoverage_MethRate_SampleCorrelation
#' 
#' @param x data for hist
#' @param ... data for hist
#'
#' @import graphics
#' 
#' @export
#'
panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="grey", ...)
}

#' Utility for HighCoverage_MethRate_SampleCorrelation
#' 
#' @param x x variable
#' @param y y variable
#' @param digits number of digits
#' @param prefix string
#' @param cex.cor graphical param
#'
#' @export
#'
#' @import graphics
#' @importFrom IRanges cor
#' @importFrom stats cor.test symnum
#' 
panel.cor <- function(x, y, digits=2, prefix="", cex.cor) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(as.numeric(cor.test(x,y)$estimate), 3)
  if(missing(cex.cor)) cex <- 0.8/strwidth(r)
  text(0.5, 0.5, r, cex=0.6*cex)
}

