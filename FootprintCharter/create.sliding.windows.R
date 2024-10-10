#' Turn sparse single molecule matrix to dense
#'
#' @import magrittr
#' @importFrom testthat equals
#' 
make.SM.mat.dense = function(binary.matrix){
  
  library(magrittr)
  
  binary.matrix %>%
    as.matrix() %>% 
    replace(equals(.,0), NA) %>% 
    replace(equals(.,1), 0) %>% 
    replace(equals(.,2), 1) -> dense.binary.matrix
  
  return(dense.binary.matrix)
  
}

#' 
#' @param single.molecule.matrix single molecule matrix (dense)
#' @param ranges.to.cover GRanges covering the loci of interest (e.g. TFBSs). Only the reads that cover continuously the entirety of the region containing these loci will be retained
#' @param params list of parameters to pass to lower functions
#' 
filter.single.molecule.matrix = function(single.molecule.matrix, ranges.to.cover=NULL, params){
  
  # Unpack parameters ####
  NA.thr.rows = params$filter.single.molecule.matrix[["NA.thr.rows"]]
  NA.thr.cols = params$filter.single.molecule.matrix[["NA.thr.cols"]]
  nanopore.pairs.filtering = params$misc[["nanopore.pairs.filtering"]]
  ########################
  
  message("### FILTERING ###")
  
  flank.width = 15
  if(nanopore.pairs.filtering){
    overall.range = IRanges::resize(ranges.to.cover, fix = "center", width = width(ranges.to.cover) + (flank.width*2))
  } else if(!is.null(ranges.to.cover)){
    overall.range = GRanges(unique(seqnames(ranges.to.cover)), IRanges(min(start(ranges.to.cover))-flank.width, max(end(ranges.to.cover))+flank.width))
  } else {
    message("no ranges.to.cover passed...proceeding with all molecules that overlap a central window")
    overall.range = GRanges("seqname", IRanges(as.integer(colnames(single.molecule.matrix)[ceiling((ncol(single.molecule.matrix)/2))]), width = 1)) %>% IRanges::resize(width = 80, fix = "center")
  }
  findOverlaps(GRanges(unique(seqnames(overall.range)), IRanges(as.integer(colnames(single.molecule.matrix)), width = 1)), overall.range) %>%
    queryHits() %>%
    unique() -> important.columns # these are cytosines we definitely want covered (e.g. TFBSs)
  reads.to.keep = rowSums(!is.na(single.molecule.matrix[,important.columns,drop=FALSE]))/length(important.columns) == 1
  single.molecule.matrix = single.molecule.matrix[reads.to.keep,,drop=FALSE]
  message(paste0("I. filtering out ", sum(!reads.to.keep), " of ", length(reads.to.keep), " reads that are discontinuous over the ranges.to.cover region"))
  
  cytosines.to.keep = colSums(is.na(single.molecule.matrix))/nrow(single.molecule.matrix) <= NA.thr.cols
  single.molecule.matrix = single.molecule.matrix[,cytosines.to.keep,drop=FALSE]
  message(paste0("II. filtering out ", sum(!cytosines.to.keep), " of ", length(cytosines.to.keep), " cytosines that are NA in more than ", NA.thr.cols*100, "% of the reads"))
  
  reads.to.keep = rowSums(is.na(single.molecule.matrix))/ncol(single.molecule.matrix) <= NA.thr.rows
  single.molecule.matrix = single.molecule.matrix[reads.to.keep,,drop=FALSE]
  message(paste0("III. filtering out ", sum(!reads.to.keep), " of ", length(reads.to.keep), " reads that have NAs in more than ", NA.thr.rows*100, "% of the remaining cytosines"))
  
  return(single.molecule.matrix)  
  
}

#' Allows to retain info about the the genomic distance between cytosines
#' 
#' @param params list of parameters to pass to lower functions
#' 
make.full.coordinate.matrix = function(single.molecule.matrix, params){
  
  # Unpack parameters ####
  padding = params$make.full.coordinate.matrix[["padding"]]
  ########################
  
  full.coordinates = (min(as.integer(colnames(single.molecule.matrix)))-padding):(max(as.integer(colnames(single.molecule.matrix)))+padding)
  missing.coordinates = full.coordinates[!full.coordinates %in% colnames(single.molecule.matrix)]
  
  cbind(
    matrix(ncol = length(missing.coordinates), nrow = nrow(single.molecule.matrix), dimnames = list(rownames(single.molecule.matrix), missing.coordinates)),
    single.molecule.matrix
  ) -> full.coordinates.matrix
  
  full.coordinates.matrix = full.coordinates.matrix[,as.character(full.coordinates),drop=FALSE]
  
  return(full.coordinates.matrix)
  
}

#' This is the core function which computes the sliding window average over single molecules to turn the data from binary to continuous
#' 
#' @param params list of parameters to pass to lower functions
matrix.sliding.window.average = function(mat, params){
  
  # Unpack parameters ####
  window.size = params$matrix.sliding.window.average[["window.size"]]
  ########################
  
  Reduce(cbind,
         mclapply(seq(ncol(mat)-window.size), function(i){
           rowMeans(mat[,i:(i+window.size-1),drop=FALSE], na.rm = TRUE)
         }, mc.cores = 1)) -> averaged.mat # on 1 core this is ~5x faster than zoo:rollmean
  
  colnames(averaged.mat) = as.integer(colnames(mat)[seq(ncol(averaged.mat))]) + ceiling(window.size/2)
  
  return(averaged.mat)
  
}

#' PROBLEM: when computing the sliding window matrix sometimes there are colummns that are etnirely NAs
#'          because two nearest cytosines are further away from each other than the size of the sliding window used.
#'          As a ducktape solution to keep working I'll fill these columns with the mean value from the two surorunding cytosines
#'          (as long as the number of columns to be completely NAs is lower than @param max.nr.cols)
#' @param sliding.window.matrix coming from a matrix.sliding.window.average call
#' @param params list of parameters to pass to lower functions
#' 
#' @importFrom dplyr first last
fill.empty.columns = function(sliding.window.matrix, params){
  
  # Unpack parameters ####
  max.nr.cols = params$fill.empty.columns[["max.nr.cols"]]
  ########################
  
  cols.to.fill = which((colSums(is.na(sliding.window.matrix))/nrow(sliding.window.matrix)) == 1)
  if (length(cols.to.fill) == 0){
    return(sliding.window.matrix)
  }
  
  split(as.integer(names(cols.to.fill)), cumsum(c(1, diff(as.integer(names(cols.to.fill))) != 1))) -> consecutive.coords

  for (i in seq_along(consecutive.coords)){
    
    current.cols.to.fill = consecutive.coords[[i]]
    
    if(dplyr::first(current.cols.to.fill) == dplyr::first(colnames(sliding.window.matrix)) | 
       dplyr::last(current.cols.to.fill) == dplyr::last(colnames(sliding.window.matrix))){
      message(paste0(length(current.cols.to.fill), " columns that need to be filles fall at the edges of the region...skipping"))
    } else if (length(current.cols.to.fill) > max.nr.cols){
      message(paste0(length(current.cols.to.fill), " columns need to be filled...more than what allowed...skipping"))
    } else {
      filling.values = rowMeans(sliding.window.matrix[,c(as.character(min(current.cols.to.fill)-1), as.character(max(current.cols.to.fill)+1)),drop=FALSE], na.rm = TRUE)
      filling.matrix = matrix(rep(filling.values, length(current.cols.to.fill)), ncol = length(current.cols.to.fill), byrow = FALSE)
      sliding.window.matrix[,as.character(current.cols.to.fill)] = filling.matrix
      message(paste0("We filled ", paste(length(consecutive.coords[[i]]), collapse = ","), " columns that were missing"))
    }
  }
  
  return(sliding.window.matrix)
  
}

#' wrapper function to go from single molecule matrix, as output by a
#' SingleMoleculeFootprinting::CallContextMethylation call, to a sliding window matrix
#' 
#' @param params list of parameters to pass to lower functions
#' 
compute.sliding.windows = function(binary.matrix, TFBSs.to.cover, params){

  binary.matrix %>%
    make.SM.mat.dense() %>%
    filter.single.molecule.matrix(., ranges.to.cover = TFBSs.to.cover, params = params) %>%
    make.full.coordinate.matrix(., params = params) %>%
    matrix.sliding.window.average(., params = params) %>%
    fill.empty.columns(., params = params) -> sliding.window.matrix.filled
  
  return(sliding.window.matrix.filled)
  
}

#' This is used to recompute sliding windows without applying the modifications necessary to counteract sparsity, which is necessary for clustering
compute.full.read.sliding.windows = function(binary.matrix, params){
  
  make.SM.mat.dense(binary.matrix) %>%
    make.full.coordinate.matrix(., params = params) %>%
    matrix.sliding.window.average(., params = params)  -> sliding.windows.matrix
  
  return(sliding.windows.matrix)
  
}


