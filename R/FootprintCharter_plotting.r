#' Plot bulk SMF separately for each partition alognside footprint detection results.
#'
#' @param MethSM sparse MethSM as returned by CallContextMethylation()
#' @param partitioned.molecules vector of partition assignments per molecule as returned by FootprintCharter()
#' @param footprints.df data.frame of footprints as returned by FootprintCharter()
#' @param TFBSs TFBSs annotation at the RegionOfInterest (optional).
#'
#' @import dplyr
#'
#' @export
#' 
#' @examples
#' 
#' Methylation = qs::qread(system.file("extdata", "Methylation_4.qs", package="SingleMoleculeFootprinting"))
#' MethSM = Methylation[[2]]
#' RegionOfInterest = GenomicRanges::GRanges("chr6", IRanges::IRanges(88106000, 88106500))
#' RegionOfInterest = IRanges::resize(RegionOfInterest, 80, "center")
#' 
#' FootprintCharter(
#'   MethSM = MethSM,
#'   RegionOfInterest = RegionOfInterest,
#'   coverage = 30,
#'   k = 16,
#'   n = 5,
#'   TF.length = c(5,75),
#'   nucleosome.length = c(120,1000),
#'   cytosine.coverage.thr = 5,
#'   verbose = TRUE
#'   ) -> FC_results
#'   
#' PlotFootprints(
#'   MethSM = MethSM, 
#'   partitioned.molecules = FC_results$partitioned.molecules, 
#'   footprints.df = FC_results$footprints.df, 
#'   TFBSs = TFBSs
#'   )
#'
PlotFootprints = function(MethSM, partitioned.molecules, footprints.df, TFBSs){
  
  partitions.id = unique(unlist(lapply(partitioned.molecules, names), use.names = FALSE, recursive = FALSE))
  lapply(partitions.id, function(i){
    lapply(partitioned.molecules, function(x){
      x[[i]]
    }) %>% unlist(use.names = FALSE)
  }) -> partitioned.molecules_list
  names(partitioned.molecules_list) = partitions.id
  
  footprints.df %<>%
    group_by(partition.nr, start, end, biological.state) %>%
    summarise(partition.coverage = sum(partition.coverage), .groups = "drop")
  
  MethSM_pooled = Reduce(rbind.fill.Matrix, MethSM)
  MethSM_dense = MethSM.to.dense(MethSM_pooled)
  lapply(partitioned.molecules_list, function(single.molecules){
    MethSM_dense[single.molecules,,drop=FALSE]
  }) -> MethSM_list
  
  if(!is.null(TFBSs)){
    TFBSs %>%
      as_tibble() %>%
      dplyr::select(start, end, TF) -> TFBS_PlottingDF
  } else {
    TFBS_PlottingDF = NULL
  }
  
  Reduce(rbind, lapply(seq_along(MethSM_list), function(i){
           data.frame(
             bulk.SMF = colSums(1 - MethSM_list[[i]], na.rm = TRUE)/colSums(!is.na(MethSM_list[[i]])) %>% as.numeric(),
             partition.nr = names(MethSM_list)[i],
             coordinate = colnames(MethSM_list[[i]]) %>% as.numeric()
           )
         })) %>%
    na.omit() %>%
    ggplot(aes(coordinate, bulk.SMF)) +
    geom_rect(data = footprints.df, mapping = aes(xmin = start-1, xmax = end+1, ymin = 0, ymax = 1, fill=biological.state), alpha=1, inherit.aes = FALSE) +
    geom_line() +
    geom_point() +
    facet_wrap(partition.nr~.) +
    geom_hline(yintercept = 0.5, linetype = 2) +
    {if(!is.null(TFBS_PlottingDF)){geom_rect(TFBS_PlottingDF, mapping = aes(xmin=start, xmax=end, ymin=-0.3, ymax=-0.1), alpha=0.5, inherit.aes = FALSE)}} +
    {if(!is.null(TFBS_PlottingDF)){geom_text(TFBS_PlottingDF, mapping = aes(x=start+((end-start)/2), y=-0.35, label=TF), inherit.aes = FALSE)}} +
    ylim(c(-0.5, 1)) +
    scale_fill_manual(breaks = c("TF", "accessible", "nucleosome", "noise", "unrecognized"), values = c("dodgerblue4", "darkseagreen", "deepskyblue3", "red3", "orange3")) +
    theme_classic()
  
}

#' Plot single molecule heatmaps of footprint detection results
#' 
#' @param footprints.df data.frame of footprints as returned by FootprintCharter()
#' @param RegionOfInterest GRanges interval to plot
#' @param partitions.order integer vector specifying the order in which to plot partitions
#' 
#' @import dplyr ggplot2
#' 
#' @export
#' 
#' @examples
#' 
#' Methylation = qs::qread(system.file("extdata", "Methylation_4.qs", package="SingleMoleculeFootprinting"))
#' MethSM = Methylation[[2]]
#' RegionOfInterest = GenomicRanges::GRanges("chr6", IRanges::IRanges(88106000, 88106500))
#' RegionOfInterest = IRanges::resize(RegionOfInterest, 80, "center")
#' 
#' FootprintCharter(
#'   MethSM = MethSM,
#'   RegionOfInterest = RegionOfInterest,
#'   coverage = 30,
#'   k = 16,
#'   n = 5,
#'   TF.length = c(5,75),
#'   nucleosome.length = c(120,1000),
#'   cytosine.coverage.thr = 5,
#'   verbose = TRUE
#'   ) -> FC_results
#'   
#' partitions.order = c(3,1,2,5,6,7,4,8)
#' Plot_FootprintCharter_SM(
#'   footprints.df = FC_results$footprints.df, 
#'   RegionOfInterest = IRanges::resize(RegionOfInterest, 500, "center"), 
#'   partitions.order = partitions.order
#'   )
#' 
Plot_FootprintCharter_SM = function(footprints.df, RegionOfInterest, partitions.order){
  
  if(is.null(partitions.order)){partitions.order = seq(max(footprints.df$partition.nr))}
  
  footprints.df %>%
    rowwise() %>%
    mutate(start = list(seq(start, end))) %>%
    ungroup() %>%
    dplyr::select(start, biological.state, partition.nr, partition.coverage, sample) %>% 
    mutate(partition.nr = factor(partition.nr, levels = rev(partitions.order))) %>%
    arrange(partition.nr) -> fp.df
  
  fp.df %>%
    dplyr::select(partition.nr, partition.coverage, sample) %>%
    distinct() %>%
    group_by(sample) %>%
    mutate(partition.coverage.cumsum = cumsum(partition.coverage), d = partition.coverage.cumsum - partition.coverage + 1) %>%
    rowwise() %>%
    mutate(read.idx = list(seq(d,partition.coverage.cumsum))) %>%
    ungroup() %>%
    dplyr::select(partition.nr, sample, read.idx) -> read.idx.df
  
  left_join(fp.df, read.idx.df, by = c("sample", "partition.nr")) %>%
    unnest(read.idx) %>%
    unnest(start) %>%
    mutate(read.idx = as.factor(read.idx)) -> PlottingDF

  PlottingDF %>%
    ggplot(aes(x = start, y = read.idx, fill=biological.state)) +
    geom_tile() +
    facet_wrap(~.data$sample, scales = "free_y", dir = 'v') +
    scale_x_continuous(limits = x.axis.breaks, breaks = x.axis.breaks, labels = format(x.axis.breaks, nsmall=1, big.mark=",")) +
    scale_fill_manual(breaks = c("TF", "accessible", "nucleosome", "noise", "unrecognized"), values = c("dodgerblue4", "darkseagreen", "deepskyblue3", "red3", "orange3")) +
    xlab(as.character(seqnames(RegionOfInterest))) + 
    ylab(paste0(sum(distinct(fp.df, partition.nr, partition.coverage)$partition.coverage), " molecules")) +
    theme_classic() +
    theme(text = element_text(size = 18), axis.title.y = element_text(vjust = 0), axis.title.x = element_text(vjust = 5),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.border = element_rect(linewidth = 1, fill = "transparent"))
  
}