#' Detects TF and nucleosome footprints enriched in a single partition
#'
#' @param MethSM sparse MethSM as returned by CallContextMethylation()
#' @param TF.length vector of two integers for footprint length bounds. Defaults to c(5,75). 
#' @param nucleosome.length vector of two integers for footprint length bounds. Defaults to c(120,1000). 
#' @param cytosine.coverage.thr Cytosine coverage threshold for footprint detection. Individual cytosines will be discarded, not whole footprints. Defaults to 5.
#' 
#' @importFrom miscTools colMedians
#' @import dplyr
#' 
.detect.footprints = function(
    SM.matrix,
    TF.length = c(5,75), 
    nucleosome.length = c(120,1000),
    cytosine.coverage.thr = 5
    ){
  
  cytosine.coverage = as.integer(apply(SM.matrix, 2, function(x){sum(!is.na(x))}))
  SM.matrix = SM.matrix[,cytosine.coverage>=cytosine.coverage.thr,drop=FALSE]
  
  median.molecule = as.numeric(1 - miscTools::colMedians(SM.matrix, na.rm = TRUE))
  median.molecule = ceiling(median.molecule) # rounds up to integer 0.5 cases
  fp.rle = rle(median.molecule)
  data.frame(
    nr.cytosines = fp.rle$lengths, 
    occupancy = fp.rle$values, 
    start = as.numeric(colnames(SM.matrix)[c(0, head(cumsum(fp.rle$lengths), -1)) + 1]),
    end = as.numeric(colnames(SM.matrix)[cumsum(fp.rle$lengths)])
    ) -> rle.df
  
  # extend footprints/accessible stretches to the middle point between cytosines
  if(nrow(rle.df)>1){
    idx.s = seq(2,nrow(rle.df))
    idx.e = seq(1,nrow(rle.df)-1)
    rle.df$new.start = rle.df$start
    rle.df$new.end = rle.df$end
    rle.df[idx.e,"new.end"] = rle.df[idx.e,"end"] + pmax(0, floor(((rle.df[idx.e+1,"start"] - rle.df[idx.e,"end"]))/2) - 1)
    rle.df[idx.s,"new.start"] = rle.df[idx.s,"start"] - pmax(0, ceiling(((rle.df[idx.s,"start"] - rle.df[idx.s-1,"end"]))/2))
    rle.df$end = rle.df$new.end
    rle.df$start = rle.df$new.start
    rle.df$new.end = NULL
    rle.df$new.start = NULL
  }
  
  rle.df %>% 
    mutate(width = end - start + 1) %>%
    mutate(biological.state = 
             ifelse(occupancy == 1 & width < TF.length[1], "noise",
             ifelse(occupancy == 1 & width >= TF.length[1] & width <= TF.length[2], "TF",
             ifelse(occupancy == 1 & width >= nucleosome.length[1] & width <= nucleosome.length[2], "nucleosome", 
             ifelse(occupancy == 0, "accessible", 
             "unrecognized"
             ))))) %>%
    mutate(biological.state = factor(biological.state, levels = c("TF", "accessible", "nucleosome", "noise", "unrecognized"))) -> footprints.df
  
  # discard TF footprints detected at the extremities of molecules 
  # where accessible cytosines are missing from one of the two flanks
  # and the true length of the footprint cannot be established
  if (footprints.df[1,"biological.state"] == "TF"){
    footprints.df[1,"biological.state"] = "unrecognized"
  }
  if (footprints.df[nrow(footprints.df),"biological.state"] == "TF"){
    footprints.df[nrow(footprints.df),"biological.state"] = "unrecognized"
  }
  
  return(footprints.df)
  
}

#' Wrapper to run the function detect.footprint across all clusters computed over a single locus
#' 
#' @param MethSM sparse MethSM as returned by CallContextMethylation()
#' @param partitioned.molecules vector of partition assignments per molecule as returned by cluster::pam()
#' @param TF.length vector of two integers for footprint length bounds. Defaults to c(5,75). 
#' @param nucleosome.length vector of two integers for footprint length bounds. Defaults to c(120,1000). 
#' @param cytosine.coverage.thr Cytosine coverage threshold for footprint detection. Individual cytosines will be discarded, not whole footprints. Defaults to 5.
#' 
#' @import dplyr
#' 
DetectFootprints = function(
    MethSM, 
    partitioned.molecules, 
    TF.length = c(5,75), 
    nucleosome.length = c(120,1000),
    cytosine.coverage.thr = 5
    ){
  
  MethSM_dense = MethSM.to.dense(MethSM)
  
  lapply(
    split(names(partitioned.molecules), partitioned.molecules), 
    function(single.molecules){MethSM_dense[single.molecules,,drop=FALSE]}
    ) -> MethSM_list
  
  Reduce(rbind, lapply(
    seq_along(MethSM_list), 
    function(i){
      .detect.footprints(
        SM.matrix = MethSM_list[[i]], 
        TF.length = TF.length, 
        nucleosome.length = nucleosome.length,
        cytosine.coverage.thr = cytosine.coverage.thr
        ) %>%
        mutate(partition.nr = names(MethSM_list)[i], partition.coverage = nrow(MethSM_list[[i]]))
         })) %>%
    arrange(partition.nr, start) %>%
    dplyr::select(partition.nr, start, end, width, partition.coverage, nr.cytosines, biological.state) -> footprints.df
  
  return(footprints.df)
  
}

#' Annotate detected TF footprints with user-provided TF motif annotations
#'
#' @param footprints.df data.frame of footprints as returned by FootprintCharter() or internally by DetectFootprints()
#' @param chromosome chromosome of current Region of interest.
#' @param TFBSs TF motif annotations. GRanges with at least two metadata columns: TF and absolute.idx for TF identity and motif index, respectively
#' 
#' @importFrom dplyr arrange mutate
#' @importFrom plyranges filter
#' @importFrom GenomicRanges GRanges findOverlaps
#'
AnnotateFootprints = function(footprints.df, chromosome, TFBSs){

  footprints.df %>%
    dplyr::arrange(start, desc(width)) %>%
    dplyr::mutate(seqnames = chromosome, TF = NA, TF.name = NA) %>%
    GRanges() -> footprints.gr
  
  TF.footprints = plyranges::filter(footprints.gr, biological.state == "TF")
  non.TF.footprints = plyranges::filter(footprints.gr, biological.state != "TF")

  # N.b.: by resizing to 1, we require that at least half the TFBSs falls into the footprint
  Overlaps = findOverlaps(TF.footprints, IRanges::resize(TFBSs, 1, "center"))
  TF.footprints$TF[unique(queryHits(Overlaps))] = split(as.character(TFBSs$TF)[subjectHits(Overlaps)], queryHits(Overlaps))
  TF.footprints$TF.name[unique(queryHits(Overlaps))] = split(as.character(TFBSs$absolute.idx)[subjectHits(Overlaps)], queryHits(Overlaps))
  footprints.df = as.data.frame(sort(c(TF.footprints, non.TF.footprints)))

  return(footprints.df)

}

#' Gather equivalent footprints by overlaps (and TF identity) under the same index
#'
#' assigns an index to footprints which allows to consider two slightly different footprints as equivalent
#' given the following condition:
#' the footprints coordinates overlap by >= 75% of the width of the smaller one and (if given) have identical TF identity
#'
#' @param footprints.df data.frame of footprints as returned by FootprintCharter() or internally by DetectFootprints() or AnnotateFootprints()
#' 
#' @importFrom dplyr  mutate arrange filter
#' @importFrom GenomicRanges GRanges elementMetadata
#' @importFrom plyranges find_overlaps
#' 
AggregateFootprints = function(footprints.df){

  groups = unique(paste(footprints.df$biological.state, footprints.df$TF, sep="."))
  biological.state.str.match = paste(paste0(unique(footprints.df$biological.state), "\\."), collapse = "|")

  lapply(seq_along(groups), function(i){

    footprints.df %>%
      mutate(TF = ifelse(is.na(TF), "NA", TF)) %>%
      arrange(desc(width)) %>%
      filter(biological.state == gsub("\\..*$", "", groups[i]) & TF == gsub(biological.state.str.match, "", groups[i], perl = TRUE)) -> footprints.df.group
    footprints.df.group$footprint.idx = NA

    GRanges(footprints.df.group) -> footprints.GRanges
    elementMetadata(footprints.GRanges) = NULL
    footprints.GRanges$original.idx = seq_along(footprints.GRanges)

    # SAFETY LOCK #
    iteration.without.progress = 0
    length.footprints.GRanges = length(footprints.GRanges)
    ###############
    while(length(footprints.GRanges) > 0){

      footprints.GRanges$flying.idx = seq_along(footprints.GRanges)
      plyranges::find_overlaps(footprints.GRanges[1], footprints.GRanges) -> overlaps
      overlaps.width = width(pintersect(footprints.GRanges[overlaps$flying.idx.x], footprints.GRanges[overlaps$flying.idx.y]))
      percentOverlap.y <- overlaps.width / width(footprints.GRanges[overlaps$flying.idx.y])
      overlaps <- overlaps[percentOverlap.y >= 0.75] # filter overlaps for a minimum percentage of overlap
      footprints.df.group[overlaps$original.idx.y, "footprint.idx"] = overlaps$original.idx.x
      footprints.GRanges[overlaps$flying.idx.y] = NULL

      # SAFETY LOCK #
      iteration.without.progress = ifelse(length(footprints.GRanges) == length.footprints.GRanges, iteration.without.progress+1, 0)
      length.footprints.GRanges = length(footprints.GRanges)
      if(iteration.without.progress > 500){stop("Stuck in infinite loop")}
      ###############

    }

    footprints.df.group %<>% arrange(start, end)

    return(footprints.df.group)

  }) -> footprints.df.group.list

  current.max.idx = 0
  for(i in seq_along(footprints.df.group.list)){

    footprints.df.group.list[[i]]$footprint.idx = footprints.df.group.list[[i]]$footprint.idx + current.max.idx
    current.max.idx = max(footprints.df.group.list[[i]]$footprint.idx)

  }

  Reduce(rbind, footprints.df.group.list) -> footprints.df

  return(footprints.df)

}