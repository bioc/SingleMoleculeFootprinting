is.binary = function(single.molecule.matrix){
  length(names(table(make.SM.mat.dense(single.molecule.matrix)))) == 2
}

#' Conditions for a Transcription Factor footprint to be identified
#'     0. footprints are detected on the "median molecule" of the passed SM.matrix 
#'     1. the footprint has to be flanked at both sides by methylated cytosines (i.e. non footprinted DNA, i.e. we can't know what lies beyond what is sequenced)
#'     2. the lengths imposed to footprints are specified in the parameters listed in params
#'
#' @param params
#' 
#' @importFrom miscTools colMedians
#' @import dplyr
#' 
detect.footprint = function(SM.matrix, params = params){
  
  # Unpack parameters ####
  TF.fp.max.length = params$detect.footprint[["TF.fp.max.length"]]
  TF.fp.min.length = params$detect.footprint[["TF.fp.min.length"]]
  Nucleosome.fp.min.length = params$detect.footprint[["Nucleosome.fp.min.length"]]
  Nucleosome.fp.max.length = params$detect.footprint[["Nucleosome.fp.max.length"]]
  cytosine.min.coverage = params$detect.footprint[["cytosine.min.coverage"]]
  ########################

  cytosine.coverage = as.integer(apply(SM.matrix, 2, function(x){sum(!is.na(x))}))
  SM.matrix = SM.matrix[,cytosine.coverage>=cytosine.min.coverage,drop=FALSE]
  
  median.molecule = as.numeric(1 - miscTools::colMedians(SM.matrix, na.rm = TRUE))
  median.molecule = ceiling(median.molecule) # with this, I'll turn cytosines that are occupied exactly half the time into fully occupied (better to be conservative)
  fp.rle = rle(median.molecule)
  data.frame(
    nr.cytosines = fp.rle$lengths, 
    occupancy = fp.rle$values, 
    start = as.numeric(colnames(SM.matrix)[c(0, head(cumsum(fp.rle$lengths), -1)) + 1]),
    end = as.numeric(colnames(SM.matrix)[cumsum(fp.rle$lengths)])
    ) -> rle.df
  # with the following 4 lines I favor occupancy stretches: i.e. when I transition from an occupied to accessible state, the occupied takes the genomic space until the next cytosine 
  # idxs = which(rle.df$occupancy == 1)
  # idxs = idxs[which(!idxs %in% c(1, nrow(rle.df)))]
  # rle.df[idxs,"start"] = rle.df[idxs-1,"end"]+1
  # rle.df[idxs,"end"] = rle.df[idxs+1,"start"]-1
  # With the following chunck I extend footprints/acc.stretches fairly, ie to the middle point between Cs
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
  
  
  #
  rle.df %>% 
    mutate(width = end - start + 1) %>%
    mutate(biological.state = 
             ifelse(occupancy == 1 & width < TF.fp.min.length, "noise",
             ifelse(occupancy == 1 & width >= TF.fp.min.length & width <= TF.fp.max.length, "TF",
             ifelse(occupancy == 1 & width >= Nucleosome.fp.min.length & width <= Nucleosome.fp.max.length, "nucleosome", 
             ifelse(occupancy == 0, "accessible", 
             "unrecognized"
             ))))) %>%
    mutate(biological.state = factor(biological.state, levels = c("TF", "accessible", "nucleosome", "noise", "unrecognized"))) -> fp.rle.df
  
  # To apply condition 1., the TF footprints detected at the extremities of the sequenced regions will be turned to "unrecognized"
  if (fp.rle.df[1,"biological.state"] == "TF"){
    fp.rle.df[1,"biological.state"] = "unrecognized"
  }
  if (fp.rle.df[nrow(fp.rle.df),"biological.state"] == "TF"){
    fp.rle.df[nrow(fp.rle.df),"biological.state"] = "unrecognized"
  }
  
  return(fp.rle.df)
  
}

#' @return single.molecule.matrix densified (if binary) and split by clustering results
#'  
prepare.footprint.detection.input = function(single.molecule.matrix, clustering.results){
  
  if (!is.binary(single.molecule.matrix)){
    warning("The single molecule matrix is NOT binary, smoothing may bias the footprint detection")
  } else {
    single.molecule.matrix %>%
      make.SM.mat.dense() -> single.molecule.matrix
  }
  
  lapply(split(names(clustering.results), clustering.results), 
         function(single.molecules){single.molecule.matrix[single.molecules,,drop=FALSE]}) -> single.molecule.matrixes_list
  
  return(single.molecule.matrixes_list)
  
}

#' Wraps detect.footprint for all clusters computed over locus
#' 
#' @param single.molecule.matrix
#' @param clustering.results
#' @param params
#' 
#' @import dplyr
#' 
detect.footprint.wrapper = function(single.molecule.matrix, clustering.results, params){
  
  prepare.footprint.detection.input(
    single.molecule.matrix = single.molecule.matrix, clustering.results = clustering.results
    ) -> single.molecule.matrixes_list
  
  Reduce(rbind,
         lapply(seq_along(single.molecule.matrixes_list), function(i){
           detect.footprint(SM.matrix = single.molecule.matrixes_list[[i]], params = params) %>%
             mutate(
               cluster.nr = names(single.molecule.matrixes_list)[i],
               cluster.coverage = nrow(single.molecule.matrixes_list[[i]])
               )
         })) %>%
    arrange(start, desc(width)) -> footprints.df
  
  return(footprints.df)
  
}

#' Arrange clusters by footprints for clean visualization of results
#' 
#' @param footprints.df as computed through detect.footprint.wrapper
#' 
arrange.clusters.by.footprints = function(footprints.df){
  
  desired.viz.order = c("TF", "TF+nucleosome", "accessible", "nucleosome")
  footprints.df %>%
    mutate(biological.state = factor(biological.state, levels = c("TF", "accessible", "nucleosome", "noise", "unrecognized"))) %>%
    group_by(cluster.nr, biological.state) %>%
    summarise(number.of.footprints = n(), .groups = "drop") %>%
    spread(biological.state, number.of.footprints, fill = 0, drop=FALSE) %>%
    dplyr::select(cluster.nr, TF, nucleosome, accessible) %>%
    mutate(prevalent.state = ifelse(TF > 0 & nucleosome > 0, "TF+nucleosome", 
                                    ifelse(TF > 0, "TF",
                                           ifelse(nucleosome > 0, "nucleosome", 
                                                  ifelse(accessible > 0, "accessible", 
                                                         ""))))) %>%
    mutate(prevalent.state = factor(prevalent.state, levels = desired.viz.order)) %>%
    arrange(prevalent.state, desc(TF)) -> arranged.df
  cluster.order = arranged.df$cluster.nr
    
  return(cluster.order)
  
}

#'
footprint.detection.diagnostic.plot = function(single.molecule.matrix, clustering.results, footprints.df, assigned.TFBSs){
  
  prepare.footprint.detection.input(
    single.molecule.matrix = single.molecule.matrix, clustering.results = clustering.results
  ) -> single.molecule.matrixes_list
  
  if(!is.null(assigned.TFBSs)){
    assigned.TFBSs %>%
      as_tibble() %>%
      dplyr::select(start, end, TF) -> TFBS_PlottingDF
  } else {
    TFBS_PlottingDF = NULL
  }
  
  
  Reduce(rbind,
         lapply(seq_along(single.molecule.matrixes_list), function(i){
           data.frame(
             Mean.molecule = colSums(1 - single.molecule.matrixes_list[[i]], na.rm = TRUE)/colSums(!is.na(single.molecule.matrixes_list[[i]])) %>% as.numeric(),
             Median.molecule = 1 - miscTools::colMedians(single.molecule.matrixes_list[[i]], na.rm = TRUE) %>% as.numeric(),
             cluster.nr = names(single.molecule.matrixes_list)[i],
             coordinate = colnames(single.molecule.matrixes_list[[i]]) %>% as.numeric()
           )
         })) %>%
    gather(mode, aggregated.molecule, -cluster.nr, -coordinate) %>%
    ggplot(aes(coordinate, aggregated.molecule)) +
    geom_line() +
    geom_point() +
    facet_grid(factor(cluster.nr, levels = suppressWarnings(arrange.clusters.by.footprints(footprints.df)))~mode) +
    geom_hline(yintercept = 0.5, linetype = 2) +
    geom_rect(data = footprints.df, mapping = aes(xmin = start-1, xmax = end+1, ymin = 0, ymax = 1, fill=biological.state), alpha=0.25, inherit.aes = FALSE) +
    {if(!is.null(TFBS_PlottingDF)){geom_rect(TFBS_PlottingDF, mapping = aes(xmin=start, xmax=end, ymin=-0.3, ymax=-0.1), alpha=0.5, inherit.aes = FALSE)}} +
    {if(!is.null(TFBS_PlottingDF)){geom_text(TFBS_PlottingDF, mapping = aes(x=start+((end-start)/2), y=-0.35, label=TF), inherit.aes = FALSE)}} +
    ggtitle("Visual assessment of footprinting detection accuracy", subtitle = "the biological.state is extrapolated from the Median.molecule") +
    theme_classic() +
    ylim(c(-0.5, 1)) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
}

# here implement reducing footprints by overlaps & TF identity
#' Gather equivalent footprints by overlaps (and TF identity) under the same index
#' 
#' effectively assigns an index to footprints which allows me later on to consider two maybe slightly different footprints as equivalent
#' given the following condition:
#' the footprints coordinates overlap by >= 75% of the width of the smaller one and (if given) have identical TF identity 
#' 
#' @param footprints.df as output by Annotate.footprints.wrapper()
#'
gather.equivalent.footprints = function(footprints.df){
  
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
      if(iteration.without.progress > 500){stop("Stuck in infinity loop")}
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


# Sub-optimal bits
# there is still some spurious TF footprints detected e.g. due to decreased coverage (see multiple NRF1 amplicon, on the right hand-side)[bait capture!!!]. 
#   These can probalby by solved by (1) coverage thr and/or (2) discarding them at the TF identity assignment step
# ambiguity in defining TF footprints that harbor TF pairs (e.g. multiple NRF1 amplicon). 
#   The solution must be to pass information across clusters at the TF identity assignment step






