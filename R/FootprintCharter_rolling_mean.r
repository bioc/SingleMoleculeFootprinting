#' Turn sparse single molecule matrix to dense
#'
#' @param MethSM sparse MethSM as returned by CallContextMethylation()
#'
#' @import magrittr
#' 
MethSM.to.dense = function(MethSM){
  
  MethSM %>%
    as.matrix() %>% 
    replace(magrittr::equals(.,0), NA) %>% 
    replace(magrittr::equals(.,1), 0) %>% 
    replace(magrittr::equals(.,2), 1) -> dense.MethSM
  
  return(dense.MethSM)
  
}

#' Filters dense matrix
#' 
#' @param MethSM sparse MethSM as returned by CallContextMethylation()
#' @param RegionOfInterest GRanges to analyse. Only the reads that cover continuously and entirely the range will be retained
#' @param verbose TRUE/FALSE
#' 
filter.dense.matrix = function(MethSM, RegionOfInterest, verbose = TRUE){
  
  NA.thr.rows = 0
  NA.thr.cols = 0
  padding = 15
  
  RegionOfInterest = IRanges::resize(RegionOfInterest, width(RegionOfInterest)+2*padding, "center")
  
  findOverlaps(GRanges(unique(seqnames(RegionOfInterest)), IRanges(as.integer(colnames(MethSM)), width = 1)), RegionOfInterest) %>%
    queryHits() %>%
    unique() -> important.columns # these are cytosines we definitely want covered (e.g. TFBSs)
  reads.to.keep = rowSums(!is.na(MethSM[,important.columns,drop=FALSE]))/length(important.columns) == 1
  MethSM = MethSM[reads.to.keep,,drop=FALSE]
  if(verbose){
    message(paste0("Discarding ", sum(!reads.to.keep), "/", length(reads.to.keep), " molecules that do not entirely cover the RegionOfInterest"))
    }
  
  cytosines.to.keep = colSums(is.na(MethSM))/nrow(MethSM) <= NA.thr.cols
  MethSM = MethSM[,cytosines.to.keep,drop=FALSE]
  if(verbose){
    message(paste0("Discarding ", sum(!cytosines.to.keep), "/", length(cytosines.to.keep), " cytosines that are NA in more than ", NA.thr.cols*100, "% of the reads"))
    }
  
  reads.to.keep = rowSums(is.na(MethSM))/ncol(MethSM) <= NA.thr.rows
  MethSM = MethSM[reads.to.keep,,drop=FALSE]
  if(verbose){
    message(paste0("Discarding ", sum(!reads.to.keep), "/", length(reads.to.keep), " molecules that have NAs in more than ", NA.thr.rows*100, "% of the remaining cytosines"))
    }
  
  return(MethSM)  
  
}

#' Computes rolling mean
#' 
#' @param MethSM single molecule matrix (dense)
#' @param window.size size of the window used to smooth molecules. Defaults to 40
#' @param padding padding size. Defaults to 20
#' 
matrix.sliding.window.average = function(MethSM, window.size = 40, padding = 20){
  
  # Full coordinates
  full.coordinates = (min(as.integer(colnames(MethSM)))-padding):(max(as.integer(colnames(MethSM)))+padding)
  missing.coordinates = full.coordinates[!full.coordinates %in% colnames(MethSM)]
  
  cbind(
    matrix(ncol = length(missing.coordinates), nrow = nrow(MethSM), dimnames = list(rownames(MethSM), missing.coordinates)),
    MethSM
  ) -> full.coordinates.matrix
  
  full.coordinates.matrix = full.coordinates.matrix[,as.character(full.coordinates),drop=FALSE]
  
  # Rolling mean
  Reduce(cbind,
         lapply(seq(ncol(full.coordinates.matrix)-window.size), function(i){
           rowMeans(full.coordinates.matrix[,i:(i+window.size-1),drop=FALSE], na.rm = TRUE)
         })) -> averaged.mat
  
  colnames(averaged.mat) = as.integer(colnames(full.coordinates.matrix)[seq(ncol(averaged.mat))]) + ceiling(window.size/2)
  
  return(averaged.mat)
  
}

#' Fills empty columns
#' 
#' when computing the sliding window matrix sometimes there are columns that are entirely NAs
#' because two nearest cytosines are further away from each other than the size of the sliding window used.
#' As a solution we fill these columns with the mean value from the two surrounding cytosines
#' (as long as the number of columns to be entirely NAs is lower than 20)
#'          
#' @param MethSM coming from a matrix.sliding.window.average call
#' @param verbose TRUE/FALSE
#' 
#' @importFrom dplyr first last
#' 
fill.empty.columns = function(MethSM, verbose = TRUE){
  
  max.nr.cols = 20
  
  cols.to.fill = which((colSums(is.na(MethSM))/nrow(MethSM)) == 1)
  if (length(cols.to.fill) == 0){return(MethSM)}
  
  split(as.integer(names(cols.to.fill)), cumsum(c(1, diff(as.integer(names(cols.to.fill))) != 1))) -> consecutive.coords

  for (i in seq_along(consecutive.coords)){
    
    current.cols.to.fill = consecutive.coords[[i]]
    
    if(dplyr::first(current.cols.to.fill) == dplyr::first(colnames(MethSM)) | 
       dplyr::last(current.cols.to.fill) == dplyr::last(colnames(MethSM))){
      if(verbose){
        message(paste0(length(current.cols.to.fill), " columns that need to be filles fall at the edges of the region...skipping"))
      }
    } else if (length(current.cols.to.fill) > max.nr.cols){
      if(verbose){
        message(paste0(length(current.cols.to.fill), " columns need to be filled...more than what allowed...skipping"))
      }
    } else {
      filling.values = rowMeans(MethSM[,c(as.character(min(current.cols.to.fill)-1), as.character(max(current.cols.to.fill)+1)),drop=FALSE], na.rm = TRUE)
      filling.matrix = matrix(rep(filling.values, length(current.cols.to.fill)), ncol = length(current.cols.to.fill), byrow = FALSE)
      MethSM[,as.character(current.cols.to.fill)] = filling.matrix
      if(verbose){
        message(paste0("We filled ", paste(length(consecutive.coords[[i]]), collapse = ","), " columns that were missing"))
      }
    }
  }
  
  return(MethSM)
  
}

#' Compute rolling mean
#' 
#' higher level wrapper
#' 
#' @param MethSM sparse MethSM as returned by CallContextMethylation()
#' @param RegionOfInterest GRanges to analyse. Only the reads that cover continuously and entirely the range will be retained
#' @param verbose TRUE/FALSE
#' 
#' @export
#' 
#' @examples
#' 
#' Methylation = qs::qread(system.file("extdata", "Methylation_4.qs", package="SingleMoleculeFootprinting"))
#' MethSM = Methylation[[2]]$SMF_MM_TKO_DE_
#' RegionOfInterest = GenomicRanges::GRanges("chr6", IRanges::IRanges(88106000, 88106500))
#' RegionOfInterest = IRanges::resize(RegionOfInterest, 80, "center")
#' 
RollingMean = function(MethSM, RegionOfInterest, verbose = TRUE){

  MethSM %>%
    MethSM.to.dense() %>%
    filter.dense.matrix(RegionOfInterest = RegionOfInterest, verbose = verbose) %>%
    matrix.sliding.window.average(window.size = 40, padding = 20) %>%
    fill.empty.columns(verbose = verbose) -> sliding.window.matrix.filled
  
  return(sliding.window.matrix.filled)
  
}
