#' Utility function to perform the dplyr full_join operation on GRanges object
#'
#' @param MethGR1 Methylation GRanges as output by the CallContextMethylation() function
#' @param MethGR2 Methylation GRanges as output by the CallContextMethylation() function
#' 
full.join.granges = function(MethGR1, MethGR2){
  
  GRanges(dplyr::full_join(as.data.frame(MethGR1), as.data.frame(MethGR2), by = c("seqnames", "start", "end", "width", "strand")))
  
}

#' Implementation performing a similar operation of plyr::rbind.fill.matrix but for sparseMatrix
#' 
#' @param x sparse matrix constructed using the function Matrix::sparseMatrix. Should have Dimnames and dims (e.g. when indexing drop=FALSE)
#' @param y sparse matrix constructed using the function Matrix::sparseMatrix. Should have Dimnames and dims (e.g. when indexing drop=FALSE)
#' 
#' @details N.b. only possible fill at the moment is 0
#' 
#' @importFrom Matrix rsparsematrix
#' 
#' @export 
#' 
rbind.fill.matrix.sparse = function(x,y){
  
  ymiss = colnames(x)[which(is.na(match(colnames(x),colnames(y))))]
  ybind = Matrix::rsparsematrix(nrow=as.double(nrow(y)),ncol=as.double(length(ymiss)),density = 0)
  colnames(ybind)<-ymiss
  
  xmiss = colnames(y)[which(is.na(match(colnames(y),colnames(x))))]
  xbind = Matrix::rsparsematrix(nrow=as.double(nrow(x)),ncol=as.double(length(xmiss)),density = 0)
  colnames(xbind) = xmiss
  
  if (ncol(xbind)>0){
    x = cbind2(x,xbind)
    x = x[,order(colnames(x)),drop=FALSE]
  }
  if(ncol(ybind)>0){
    y = cbind2(y,ybind)
    y = y[,order(colnames(y)),drop=FALSE]
  }
  
  result = rbind2(x,y[,order(match(colnames(y),colnames(x))),drop=FALSE])
  if (all(result@Dim > 0)){
    rownames(result) = c(rownames(x), rownames(y)) # for some reason rbind2 drop rownames
    # result = result[,order(colnames(result), decreasing = FALSE)] # This shouldn't be necessary for rows
  }
  
  return(result)
  
}

#' Implementation performing a similar operation of rbind.fill.matrix.sparse but for columns
#' 
#' @param x sparse matrix constructed using the function Matrix::sparseMatrix. Should have Dimnames and dims (e.g. when indexing drop=FALSE)
#' @param y sparse matrix constructed using the function Matrix::sparseMatrix. Should have Dimnames and dims (e.g. when indexing drop=FALSE)
#' 
#' @details N.b. only possible fill at the moment is 0
#' 
#' @export 
#' 
cbind.fill.matrix.sparse = function(x,y){
  
  ymiss = rownames(x)[which(is.na(match(rownames(x),rownames(y))))]
  ybind = rsparsematrix(nrow=as.double(length(ymiss)),ncol=as.double(ncol(y)),density = 0)
  rownames(ybind) = ymiss
  
  xmiss = rownames(y)[which(is.na(match(rownames(y),rownames(x))))]
  xbind = rsparsematrix(nrow=as.double(length(xmiss)),ncol=as.double(ncol(x)),density = 0)
  rownames(xbind) = xmiss
  
  x = rbind2(x,xbind)
  y = rbind2(y,ybind)
  
  result = cbind2(x,y[order(match(rownames(y),rownames(x))),,drop=FALSE])
  if (all(result@Dim > 0)){
    colnames(result) = c(colnames(x), colnames(y)) # for some reason cbind2 drop colnames
    result = result[,order(colnames(result), decreasing = FALSE),drop=FALSE]
  }
  
  return(result)
  
}

#' Utility function to remove cytosines whose MTase target genomic context is affected by SNPs
#' 
#' @param Methylation as output by the CallContextMethylation() function
#' @param CytosinesToMask 
#' @param MaskSMmat whether the parameter Methylation includes single molecule matrixes
#' @param SampleStringMatch list of per-sample string matches that are used to uniquely identify the relevant column for each species in the Methylation object. Defaults to list(Cast = "_CTKO", Spret = "_STKO")
#' @param Experiment as detected by the DetectExperimentType() function. Should be either "DE" or "NO"
#' 
#' @export
#'
MaskSNPs = function(Methylation, CytosinesToMask, MaskSMmat = FALSE, SampleStringMatch = list(Cast = "_CTKO", Spret = "_STKO"), Experiment){
  
  Castaneus.string.match = SampleStringMatch$Cast
  Spretus.string.match = SampleStringMatch$Spret
  
  if(Experiment == "DE"){
    message("Masking GRanges in DE mode")
    if (MaskSMmat){
      Cast_DE_disrupted = unique(queryHits(findOverlaps(Methylation[[1]], CytosinesToMask[CytosinesToMask$DisruptedInCast])))
      Spret_DE_disrupted = unique(queryHits(findOverlaps(Methylation[[1]], CytosinesToMask[CytosinesToMask$DisruptedInSpret])))
      CastCols = grep(Castaneus.string.match, colnames(elementMetadata(Methylation[[1]])))
      SpretCols = grep(Spretus.string.match, colnames(elementMetadata(Methylation[[1]])))
      if(length(CastCols)>0){elementMetadata(Methylation[[1]])[Cast_DE_disrupted,CastCols] = NA}
      if(length(SpretCols)>0){elementMetadata(Methylation[[1]])[Spret_DE_disrupted,SpretCols] = NA}
      
      message("Masking SM matrix")
      CastDisruptedCoords = start(Methylation[[1]])[Cast_DE_disrupted]
      SpretDisruptedCoords = start(Methylation[[1]])[Spret_DE_disrupted]
      CastMats = grep(Castaneus.string.match, grep("_MethRate$", colnames(elementMetadata(Methylation[[1]])), value = TRUE))
      SpretMats = grep(Spretus.string.match, grep("_MethRate$", colnames(elementMetadata(Methylation[[1]])), value = TRUE))
      
      Methylation[[2]][CastMats] = lapply(Methylation[[2]][CastMats], function(mat){mat[,!colnames(mat) %in% CastDisruptedCoords, drop=FALSE]})
      Methylation[[2]][SpretMats] = lapply(Methylation[[2]][SpretMats], function(mat){mat[,!colnames(mat) %in% SpretDisruptedCoords, drop=FALSE]})
      
    } else {
      message("Skipping SM matrix")
      Cast_DE_disrupted = unique(queryHits(findOverlaps(Methylation, CytosinesToMask[CytosinesToMask$DisruptedInCast])))
      Spret_DE_disrupted = unique(queryHits(findOverlaps(Methylation, CytosinesToMask[CytosinesToMask$DisruptedInSpret])))
      CastCols = grep(Castaneus.string.match, colnames(elementMetadata(Methylation)))
      SpretCols = grep(Spretus.string.match, colnames(elementMetadata(Methylation)))
      if(length(CastCols)>0){elementMetadata(Methylation)[Cast_DE_disrupted,CastCols] = NA}
      if(length(SpretCols)>0){elementMetadata(Methylation)[Spret_DE_disrupted,SpretCols] = NA}
    }
    
  } else if (Experiment == "NO") {
    message("Masking GRanges in NO mode")
    for(context in seq(2)){
      if (MaskSMmat){
        Cast_SE_disrupted = unique(queryHits(findOverlaps(Methylation[[1]][[context]], CytosinesToMask[CytosinesToMask$DisruptedInCast])))
        Spret_SE_disrupted = unique(queryHits(findOverlaps(Methylation[[1]][[context]], CytosinesToMask[CytosinesToMask$DisruptedInSpret])))
        CastCols = grep(Castaneus.string.match, colnames(elementMetadata(Methylation[[1]][[context]])))
        SpretCols = grep(Spretus.string.match, colnames(elementMetadata(Methylation[[1]][[context]])))
        if(length(CastCols)>0){elementMetadata(Methylation[[1]][[context]])[Cast_SE_disrupted,CastCols] = NA}
        if(length(SpretCols)>0){elementMetadata(Methylation[[1]][[context]])[Spret_SE_disrupted,SpretCols] = NA}
        
        message("Masking SM matrix")
        CastDisruptedCoords = start(Methylation[[1]][[context]])[Cast_SE_disrupted]
        SpretDisruptedCoords = start(Methylation[[1]][[context]])[Spret_SE_disrupted]
        CastMats = grep(Castaneus.string.match, grep("_MethRate$", colnames(elementMetadata(Methylation[[1]][[context]])), value = TRUE))
        SpretMats = grep(Spretus.string.match, grep("_MethRate$", colnames(elementMetadata(Methylation[[1]][[context]])), value = TRUE))
        
        for(sample in CastMats){
          Coordinates_to_keep = !colnames(Methylation[[2]][[sample]][[context]]) %in% CastDisruptedCoords
          Methylation[[2]][[sample]][[context]] = Methylation[[2]][[sample]][[context]][,Coordinates_to_keep, drop=FALSE]
        }
        
        for(sample in SpretMats){
          Coordinates_to_keep = !colnames(Methylation[[2]][[sample]][[context]]) %in% SpretDisruptedCoords
          Methylation[[2]][[sample]][[context]] = Methylation[[2]][[sample]][[context]][,Coordinates_to_keep, drop=FALSE]
        }
      } else {
        message("Skipping SM matrix")
        Cast_SE_disrupted = unique(queryHits(findOverlaps(Methylation[[context]], CytosinesToMask[CytosinesToMask$DisruptedInCast])))
        Spret_SE_disrupted = unique(queryHits(findOverlaps(Methylation[[context]], CytosinesToMask[CytosinesToMask$DisruptedInSpret])))
        CastCols = grep(Castaneus.string.match, colnames(elementMetadata(Methylation[[context]])))
        SpretCols = grep(Spretus.string.match, colnames(elementMetadata(Methylation[[context]])))
        if(length(CastCols)>0){elementMetadata(Methylation[[context]])[Cast_SE_disrupted,CastCols] = NA}
        if(length(SpretCols)>0){elementMetadata(Methylation[[context]])[Spret_SE_disrupted,SpretCols] = NA}
      }
    }
  } else {stop("Unrecognized Experiment")}
  
  return(Methylation)
  
}

#' Convenience function to arrange a list of given TFBSs into clusters
#' 
#' For each TFBS, the genomic neighborhood defined by max_cluster_width will be scanned for adjacent TFBSs. 
#' The hits will be filtered for min_intersite_distance where, in case of overlapping TFBSs, the second TFBS will be arbitrarily dropped.
#' These TFBSs plus the central "anchoring" one will define a TFBS cluster. 
#' This approach implies that the same TFBS can be employed to design multiple clusters in a sliding-window fashion.
#' 
#' @param TFBSs GRanges object of TFBSs
#' @param max_intersite_distance maximum allowed distance in base pairs between two TFBS centers for them to be considered part of the same cluster. Defaults to 75.
#' @param min_intersite_distance minimum allowed distance in base pairs between two TFBS centers for them not to be discarded as overlapping. 
#'                               This parameter should be set according to the width of the bins used for later sorting. Defaults to 15.
#' @param max_cluster_size maximum number of TFBSs to be contained in any given cluster. Defaults to 6
#' @param max_cluster_width maximum width of TFBS clusters in bps. Defaults to 300
#' @param add.single.TFs Whether to add the TFs not used to create TFBS.clusters to the list for sorting. Defaults to TRUE
#' 
#' @import GenomicRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom IRanges IRanges reduce
#' @importFrom plyranges filter_by_non_overlaps
#' 
#' @return list with two elements: ClusterCoordinates (GRanges object of clusters coordinates) and ClusterComposition (GRangesList of sites for each cluster)
#' 
#' @export
#'
Arrange_TFBSs_clusters = function(
    TFBSs, 
    max_intersite_distance = 75, min_intersite_distance = 15, 
    max_cluster_size = 6, max_cluster_width = 300, 
    add.single.TFs = TRUE){
  
  collection_window_width = max_intersite_distance*2
  TFBSs_resized_1 = resize(TFBSs,1,fix='center')
  TFBSs_resized_window = resize(TFBSs,collection_window_width,fix='center')
  
  Overlaps = findOverlaps(TFBSs_resized_window, TFBSs_resized_1, ignore.strand=TRUE)
  pair_dist = start(TFBSs_resized_1[subjectHits(Overlaps)]) - start(TFBSs_resized_1[queryHits(Overlaps)])
  
  message("Removing self-overlaps and redundant pairs")
  message("Removing pairs containing overlapping factors")
  Overlaps_filtered = Overlaps[pair_dist>min_intersite_distance,]
  
  message("Constructing GRanges object of clusters coordinates")
  TF_cluster = IRanges::reduce(GRanges(
    seqnames(TFBSs[queryHits(Overlaps_filtered)]),
    IRanges(start(TFBSs[queryHits(Overlaps_filtered)]), 
            end(TFBSs[subjectHits(Overlaps_filtered)]))
  )) # very rarely these exceed width of 300bp
  TF_cluster = sort(TF_cluster)
  
  message("Computing number of sites per cluster")
  TF_cluster$number_of_TF = countOverlaps(TF_cluster,TFBSs)
  
  message(paste0("Discaring clusters with more than ", max_cluster_size, "sites"))
  TF_cluster = TF_cluster[TF_cluster$number_of_TF <= max_cluster_size]
  
  message("Creating TFBS_cluster ID on the fly")
  names(TF_cluster) = paste0('TFBS_cluster_',seq_along(TF_cluster))
  
  message("Constructing GRangesList of sites per cluster")
  Overlaps_clusters = findOverlaps(TF_cluster,TFBSs_resized_1)
  TF_list = split(TFBSs[subjectHits(Overlaps_clusters)], queryHits(Overlaps_clusters))
  names(TF_list) = names(TF_cluster)
  
  return(list(ClusterCoordinates = TF_cluster, ClusterComposition = TF_list))
  
}

#' Create methylation calling windows to call context methylation in one run for clusters lying proximally to each other
#' 
#' Relevant for genome-wide analyses
#' 
#' @param TFBS_cluster_coordinates TFBS cluster coordinates analogous to ClusterCoordinates object returned by Arrange_TFBSs_clusters function
#' @param max_intercluster_distance maximum distance between two consecutive TFBS clusters for them to be grouped in the same window
#' @param max_window_width upper limit to window width. This value should be adjusted according to the user's system as it determines the amount of memory used in the later context methylation call
#' @param min_cluster_width lower limit to window width. Corresponds to the scenario when a window contains a single TFBS cluster.
#' @param genomic.seqlenghts used to fix the windows spanning over chromosome edges. To be fetched by GenomeInfoDb::seqlengths() or equivalent.
#' @param fix.window.size Defaults to FALSE. When TRUE, overrides arguments max_intercluster_distance and max_window_width and produces windows containing a fixed number of TFBS_clusters.
#' @param max.window.size Max number of TFBS_clusters per window. Used only when fix.window.size is TRUE. N.b.: window size could be slightly higher than passed value if TFBS_cluster_coordinates overlap
#' 
#' @import GenomicRanges
#' @import plyranges
#' 
#' @return GRanges object of window coordinates to be used for more efficient calls of CallContextMethylation 
#' 
#' @export
#' 
Create_MethylationCallingWindows = function(
    TFBS_cluster_coordinates, 
    max_intercluster_distance = 100000, 
    max_window_width = 5000000, 
    min_cluster_width = 600, 
    genomic.seqlenghts, 
    fix.window.size = FALSE, 
    max.window.size = 500){
  
  if(isFALSE(fix.window.size)){
    
    message(paste0("Group TFBS_clusters that fall within ", max_intercluster_distance, "bp from each other in broader searching windows"))
    SearchingWindows = plyranges::reduce_ranges(resize(TFBS_cluster_coordinates, width = max_intercluster_distance, fix = 'center'))
    message("Trimming searching windows")
    start(SearchingWindows) = start(SearchingWindows) + ((max_intercluster_distance/2) - (min_cluster_width/2))
    end(SearchingWindows) = end(SearchingWindows) - ((max_intercluster_distance/2) - (min_cluster_width/2))
    
    TooLarge = width(SearchingWindows) > max_window_width
    
    while(sum(TooLarge) > 0){
      
      SmallerWindow = max_intercluster_distance*0.9
      message("Reducing too large searching windows")
      
      Overlaps = findOverlaps(SearchingWindows[TooLarge], TFBS_cluster_coordinates)
      SmallerWindows = plyranges::reduce_ranges(resize(TFBS_cluster_coordinates[subjectHits(Overlaps)], width = SmallerWindow, fix = "center"))
      SearchingWindows = sort(c(SearchingWindows[!TooLarge], SmallerWindows))
      
      max_intercluster_distance = SmallerWindow
      TooLarge = width(SearchingWindows) > max_window_width
      
      message("fix windows spanning over chromosome edges")
      start(SearchingWindows)[which(start(SearchingWindows) < 0)] = 1
      for(chr in as.character(sort(unique(seqnames(SearchingWindows))))){
        end(SearchingWindows[seqnames(SearchingWindows) == chr & end(SearchingWindows) > sort(genomic.seqlenghts)[chr]]) = sort(genomic.seqlenghts)[chr]
      }
      
    }
    
  } else {
    
    TFBS_cluster_coordinates %>%
      IRanges::reduce(ignore.strand = TRUE) %>%
      sort() %>%
      plyranges::mutate(idx = rep(seq(ceiling(length(.)/max.window.size)), each = max.window.size)[seq_along(.)]) %>%
      plyranges::group_by(seqnames, idx) %>%
      plyranges::summarise(start = min(start)-1, end = max(end)+1) %>% # padding for strand collapsing 
      #       N.b. when working with overlapping GenomicTiles, 
      #            the resulting windows will be precise from the start of the first tile to the end of the last, 
      #            but will inherit the same overlap as the GenomicTiles
      GRanges() %>%
      plyranges::select(-idx) -> SearchingWindows
    
  }
  
  return(sort(SearchingWindows))
  
}
