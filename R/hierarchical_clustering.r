#' Perform Hierarchical clustering on single reads
#'
#' @param MethSM Single molecule methylation matrix
#'
#' @importFrom stats hclust
#' @importFrom stats dist
#'
HierarchicalClustering = function(MethSM){
  
  SubsetSize = 500
  if(nrow(MethSM) > SubsetSize){ #subset to 500 molecules to avoid problem with Hc
    MethSM_subset = MethSM[sample(dimnames(MethSM)[[1]],SubsetSize),]
  }else{
    SubsetSize = nrow(MethSM)
    MethSM_subset = MethSM
  }
  ReadsDist = dist(MethSM_subset)
  iteration = 0
  while(sum(is.na(ReadsDist)) > 0){ # sometimes dist between some pair of reads is NA, possibly because of no overlapping Cs
    iteration = iteration + 1
    if (iteration > SubsetSize) {
      SubsetSize = ceiling(SubsetSize*0.9)
      iteration = 0
    }
    MethSM_subset = MethSM[sample(dimnames(MethSM)[[1]],SubsetSize),]
    ReadsDist = dist(MethSM_subset)
  }
  hc=hclust(ReadsDist)
  MethSM_HC = MethSM_subset[hc$order,]
  
  return(MethSM_HC)
  
}