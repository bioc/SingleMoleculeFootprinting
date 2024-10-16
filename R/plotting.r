#' Plot average methylation
#'
#' @param MethGR Average methylation GRanges obj
#' @param MethSM Single molecule matrix(es)
#' @param RegionOfInterest GRanges interval to plot
#' @param SortedReads List of sorted reads, needs to be passed along with the parameter MethSM. If both are passed, only counts relevant to sorting will be plotted
#' @param ShowContext TRUE or FALSE (default). Causes the genomic context of the plotted cytosines to be displayed as the dot shape
#' @param TFBSs GRanges object of transcription factor binding sites to include in the plot. Assumed to be already subset. Also assumed that the tf names are under the column "TF"
#' @param SNPs GRanges object of SNPs to visualize. Assumed to be already subset. Assumed to have the reference and alternative sequences respectively under the columns "R" and "A"
#' @param SortingBins GRanges object of sorting bins (absolute) coordinate to visualize
#'
#' @import GenomicRanges
#' @import tidyverse
#' @importFrom plyr .
#' @importFrom stats na.omit
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
#'
#' @examples
#'
#' library(GenomicRanges)
#'
#' RegionOfInterest = GRanges("chr12", IRanges(20464551, 20465050))
#' Methylation = qs::qread(system.file("extdata", "Methylation_3.qs", 
#' package="SingleMoleculeFootprinting"))
#' TFBSs = qs::qread(system.file("extdata", "TFBSs_3.qs", package="SingleMoleculeFootprinting"))
#'
#' PlotAvgSMF(MethGR = Methylation[[1]], RegionOfInterest = RegionOfInterest, TFBSs = TFBSs)
#'
PlotAvgSMF = function(MethGR, MethSM=NULL, RegionOfInterest, SortedReads=NULL, ShowContext=FALSE, TFBSs=NULL, SNPs=NULL, SortingBins=NULL){

  # Prepare SMF data
  MethGR %>%
    as_tibble() %>%
    dplyr::select(-grep("_Coverage$", colnames(.)), -.data$end, -.data$width, -.data$strand) %>%
    gather(sample, MethRate, -.data$seqnames, -.data$start, -.data$GenomicContext) %>%
    na.omit() -> PlottingDF
  
  if (any(stringr::str_detect(colnames(values(MethGR)), "^A_|^R_"))){
    PlottingDF$sample = factor(PlottingDF$sample, levels = sort(unique(PlottingDF$sample), decreasing = TRUE))
  }
  
  if(!is.null(SortedReads) & !is.null(MethSM)){
    message("Sorted reads passed along with SM matrix...plotting only counts relevant to sorting")
    if(length(MethSM) != length(SortedReads)){
      stop("Number of samples in MethSM and SortedReads do not correspond...quitting")
    }
    # here I recalculate MethRate and replace it in the PlottingDF
    # removing cytosines not covered by relevant reads
    # No need to take care of coverage again since we only care for counts over sorting bins,
    # and these reads are by definition covering bins enough
    Reduce(rbind,
    lapply(seq_along(SortedReads), function(i){
      RecalculatedMethRate = colMeans_drop0(MethSM[[i]][as.character(unlist(SortedReads[[i]])),]) - 1
      PlottingDF %>%
        filter(sample == paste0(names(SortedReads)[i], "_MethRate")) %>%
        filter(start %in% names(RecalculatedMethRate)) %>%
        arrange(start) %>%
        mutate(MethRate = as.double(RecalculatedMethRate[order(names(RecalculatedMethRate))]))
    })) -> PlottingDF
  } else if (!is.null(SortedReads) & is.null(MethSM)){
    stop("Sorted reads passed without SM matrix...please provide SM matrix")
  } else if (is.null(SortedReads)){
    message("No sorted reads passed...plotting counts from all reads")
  }
  
  OurFavouriteColors = c("Black", RColorBrewer::brewer.pal(n = 9, name = "Set1"))
  ColorsToUse = OurFavouriteColors[seq_along(unique(PlottingDF$sample))]

  # Prepare TFBS
  if(!is.null(TFBSs)){
    TFBSs %>%
      as_tibble() %>%
      dplyr::select(.data$start, .data$end, .data$TF) -> TFBS_PlottingDF
  }

  # Prepare SNPs
  if(!is.null(SNPs)){
    SNPs %>%
      as_tibble() %>%
      rowwise() %>%
      mutate(width = ifelse(max(nchar(.data$R), nchar(.data$A)) == 1, 3, max(nchar(.data$R), nchar(.data$A)))) %>%
      dplyr::select(.data$start, .data$width) %>%
      mutate(y_coord = -0.13) -> SNPs_PlottingDF
  }

  # Prepare SortingBins
  if(!is.null(SortingBins)){
    SortingBins %>%
      as.data.frame() %>%
      select(.data$start, .data$end) -> Bins_PlottingDF
  }

  PlottingDF %>%
    ggplot(aes(x=.data$start, y=1-.data$MethRate, color=.data$sample)) +
    geom_line() +
    {if(ShowContext){geom_point(aes(shape=.data$GenomicContext))}else{geom_point()}} +
    {if(!is.null(TFBSs)){geom_rect(TFBS_PlottingDF, mapping = aes(xmin=.data$start, xmax=.data$end, ymin=-0.09, ymax=-0.04), inherit.aes = FALSE)}} +
    {if(!is.null(TFBSs)){ggrepel::geom_text_repel(TFBS_PlottingDF, mapping = aes(x=.data$start+((.data$end-.data$start)/2), y=-0.02, label=.data$TF), min.segment.length = .1, max.overlaps = 1e+05, inherit.aes = FALSE)}} +
    {if(!is.null(SNPs)){geom_tile(SNPs_PlottingDF, mapping = aes(x=.data$start, y=.data$y_coord, width=.data$width), color = ColorsToUse[2], fill = ColorsToUse[2], height = 0.05, inherit.aes = FALSE)}} +
    {if(!is.null(SortingBins)){geom_rect(Bins_PlottingDF, mapping = aes(xmin=.data$start, xmax=.data$end, ymin=-0.02, ymax=0), color="black", fill="white", inherit.aes = FALSE)}} +
    geom_hline(aes(yintercept=0)) +
    ylab("SMF") +
    xlab("") +
    ylim(c(-0.25,1)) +
    xlim(c(start(RegionOfInterest),end(RegionOfInterest))) +
    ggtitle(RegionOfInterest) +
    scale_color_manual(values = ColorsToUse) +
    theme_classic()

}

#' Plot single molecule stack
#'
#' @param MethSM Single molecule methylation matrix
#' @param RegionOfInterest GRanges interval to plot
#'
#' @import GenomicRanges
#' @import tidyverse
#' @importFrom tibble rownames_to_column
#' @importFrom stats na.omit
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' 
#' library(GenomicRanges)
#' 
#' RegionOfInterest = GRanges("chr12", IRanges(20464551, 20465050))
#' Methylation = qs::qread(system.file("extdata", "Methylation_3.qs", 
#' package="SingleMoleculeFootprinting"))
#'
#' PlotSingleMoleculeStack(MethSM = Methylation[[2]], RegionOfInterest = RegionOfInterest)
#'
PlotSingleMoleculeStack = function(MethSM, RegionOfInterest){
  
  Reduce(rbind, lapply(seq_along(MethSM), function(i){
    MethSM[[i]] %>%
      as.matrix() %>%
      as.data.frame() %>%
      rownames_to_column(var = "ReadID") %>%
      mutate(Sample = names(MethSM)[i]) %>%
      gather(Coordinate, Methylation, -.data$ReadID, -.data$Sample) %>%
      mutate(Methylation = ifelse(.data$Methylation == 0, NA, .data$Methylation-1))
  })) %>%
    na.omit() %>%
    mutate(Methylation = as.factor(.data$Methylation), Coordinate = as.numeric(.data$Coordinate)) -> PlottingDF
  PlottingDF$ReadID = factor(PlottingDF$ReadID, levels = Reduce(c, lapply(MethSM, rownames)))
  
  PlottingDF %>%
    group_by(.data$Sample) %>%
    summarise(NrReads = length(unique(.data$ReadID))) %>%
    ungroup() %>%
    mutate(Label = paste0(Sample, " (", .data$NrReads, " reads)")) %>%
    dplyr::select(.data$Sample, .data$Label) -> Labels
  names(Labels$Label) = Labels$Sample

  PlottingDF %>%
    ggplot(aes(x=.data$Coordinate, y=.data$ReadID)) + 
    geom_tile(aes(fill=.data$Methylation), height=1, width=5) +
    facet_wrap(~.data$Sample, scales = "free_y", dir = 'v', 
               labeller = as_labeller(Labels$Label)) +
    ylab("") +
    xlab("") +
    xlim(c(start(RegionOfInterest),end(RegionOfInterest))) +
    scale_discrete_manual(aesthetics = "fill", values = c("black", "grey")) +
    theme_classic() +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())

}

.arrange.MethSM.by.SortedReads = function(MethSM, SortedReads, ordered.sorting.patterns){
  
  NAMES = names(MethSM)
  if(!is.null(ordered.sorting.patterns)){
    MethSM = lapply(seq_along(MethSM), function(i){
      MethSM[[i]][unlist(SortedReads[[i]][ordered.sorting.patterns]),]
    })
  } else { # this because indexing with NULL (when sorting.strategy == "custom") return character(0)
    MethSM = lapply(seq_along(MethSM), function(i){
      MethSM[[i]][unlist(SortedReads[[i]]),]
    })
  }

  names(MethSM) = NAMES
  
  return(MethSM)
  
}

#' Wrapper for PlotSingleMoleculeStack function
#'
#' adds the convenience of arranging reads before plotting
#'
#' @param MethSM Single molecule methylation matrix
#' @param RegionOfInterest GRanges interval to plot
#' @param sorting.strategy One of "classical" (default), "custom", "hierarchical.clustering" or "None".
#' Set to "classical" for classical one-TF/TF-pair sorting (as described in Sönmezer et al, MolCell, 2021). Should be passed along with argument SortedReads set to the Sorted reads object as returned by SortReads function.
#' If set to "custom", SortedReads should be a list with one item per sample (corresponding to MethSM).
#' If set to "hierarchical.clustering", the function will perform hierarchical clustering in place on a subset of reads. Useful to check for duplicated reads in amplicon sequencing experiments.
#' If set to "None", it will plot unsorted reads.
#' The argument sorting,strategy will always determine how to display reads with priority over the argument SortedReads
#' @param SortedReads Defaults to NULL, in which case will plot unsorted reads. Sorted reads object as returned by SortReads function 
#'  
#'
#' @export
#'
#' @examples
#' 
#' library(GenomicRanges)
#' 
#' RegionOfInterest = GRanges("chr12", IRanges(20464551, 20465050))
#' Methylation = qs::qread(system.file("extdata", "Methylation_3.qs", 
#' package="SingleMoleculeFootprinting"))
#' TFBSs = qs::qread(system.file("extdata", "TFBSs_3.qs", 
#' package="SingleMoleculeFootprinting"))
#' SortedReads = SortReadsBySingleTF(MethSM = Methylation[[2]], TFBS = TFBSs)
#'
#' PlotSM(MethSM = Methylation[[2]], RegionOfInterest = RegionOfInterest, SortedReads = SortedReads)
#'
PlotSM = function(MethSM, RegionOfInterest, sorting.strategy="classical", SortedReads = NULL){
  
  #### 1.
  if(sorting.strategy == "classical" & all(unlist(lapply(SortedReads, is.list)))){
    
    if (length(MethSM) != length(SortedReads)){stop("Number of samples is not consistent between SortedReads and MethSM...quitting")}
    message("Arranging reads according to classical sorting.strategy")
    PatternLength = unique(unlist(lapply(seq_along(SortedReads), function(i){unique(nchar(names(SortedReads[[i]])))})))
    
    if (PatternLength == 3){ # Single TF
      message("Inferring sorting was performed by single TF")
      ordered.sorting.patterns = rev(Reduce(c, SingleTFStates()))
    } else if (PatternLength == 4){ # TF pair
      message("Inferring sorting was performed by TF pair")
      # ordered.sorting.patterns = rev(as.character(unlist(TFPairStates())))
      ordered.sorting.patterns = rev(c("1001", "1010", "1011", "0101", "1101", "1111", "0111", "0011", "0001", "1110", "0110", "0010", "1100", "0100", "1000", "0000"))
    } else {
      ordered.sorting.patterns = NULL
    }
    
    MethSM = .arrange.MethSM.by.SortedReads(MethSM, SortedReads, ordered.sorting.patterns)
    
  #### 2.
  } else if (sorting.strategy == "hierarchical.clustering"){
    
    if(!is.null(SortedReads)){warning("Ignoring passed SortedReads and performing hierarchical clustering")}
    message("Perfoming hierarchical clustering on single molecules before plotting")
    MethSM = lapply(MethSM, HierarchicalClustering)
    
  #### 3.
  } else if (sorting.strategy == "custom"){
    
    if (length(MethSM) != length(SortedReads)){stop("Number of samples is not consistent between SortedReads and MethSM...quitting")}
    message("Arranging reads according to custom sorting.strategy")
    MethSM = .arrange.MethSM.by.SortedReads(MethSM, SortedReads, ordered.sorting.patterns=NULL)
    
  #### 4.
  } else if (sorting.strategy == "None"){
    
    if(!is.null(SortedReads)){warning("Ignoring passed SortedReads and plotting unsorted reads")}
    message("No sorting passed or specified, will plot unsorted reads")
    
  #### 5.
  } else {stop("Invalid value for sorting.strategy")}
  
  PlotSingleMoleculeStack(MethSM, RegionOfInterest)

}

#' Plot states quantification bar
#'
#' @param SortedReads Sorted reads object as returned by SortReads function
#' @param states either SingleTFStates() or TFPairStates()
#' 
#' @return Bar plot quantifying states
#' 
#' @import dplyr
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' 
#' library(GenomicRanges)
#' 
#' RegionOfInterest = GRanges("chr12", IRanges(20464551, 20465050))
#' Methylation = qs::qread(system.file("extdata", "Methylation_3.qs", 
#' package="SingleMoleculeFootprinting"))
#' TFBSs = qs::qread(system.file("extdata", "TFBSs_3.qs", package="SingleMoleculeFootprinting"))
#' SortedReads = SortReadsBySingleTF(MethSM = Methylation[[2]], TFBS = TFBSs)
#'
#' StateQuantificationPlot(SortedReads = SortedReads, states = SingleTFStates())
#'
StateQuantificationPlot = function(SortedReads, states){
  
  if(states[[1]] == "101"){
    states_df = unique(gather(as_tibble(states), "State", "Pattern")) %>%
      mutate(State = factor(.data$State, levels = c("bound", "accessible", "closed", "unassigned")))
  } else if (states[[1]] == "1001"){
    states_df = data.frame(
      State = rep(names(states), lengths(states)),
      Pattern = unlist(states, use.names = FALSE)
    ) %>%
      mutate(
        State = factor(State, levels = c("bound_bound", "bound_misc", "misc_bound", "misc")),
        Pattern = factor(Pattern, levels = c("1001", "1010", "1011", "0101", "1101", "1111", "0111", "0011", "0001", "1110", "0110", "0010", "1100", "0100", "1000", "0000"))
      ) %>%
      arrange(State, Pattern)
  }
  
  Reduce(rbind, lapply(seq_along(SortedReads), function(i){
    data.frame(
      ReadID = unlist(SortedReads[[i]], recursive = TRUE, use.names = FALSE),
      Pattern = rep(names(SortedReads[[i]]), lengths(SortedReads[[i]])),
      Sample = names(SortedReads)[i]
    )
  })) %>%
    left_join(., states_df, by = "Pattern") %>%
    arrange(desc(.data$State)) %>%
    mutate(ReadID = factor(.data$ReadID, levels = .data$ReadID)) -> PlottingDF
  
  if(states[[1]] == "101"){
    PlottingDF %>%
      ggplot(aes(x=1, y=.data$ReadID, fill=.data$State)) + 
      geom_tile() +
      facet_wrap(~.data$Sample, scales = "free_y", dir = 'v') +
      ylab("") +
      xlab("") +
      scale_fill_manual(breaks = names(states), values = c("#984EA3", "#4DAF4A", "#377EB8", "#999999")) + 
      theme_classic() +
      theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.line = element_blank()) -> pl
  } else if (states[[1]] == "1001") {
    PlottingDF %>%
      separate(.data$Pattern, into = c(paste0("Bin", seq(unique(nchar(unlist(states)))))), sep = "(?<=.)", extra = 'drop') %>%
      gather(Bin, Methylation, -.data$ReadID, -.data$Sample, -.data$State) %>%
      ggplot(aes(x=.data$Bin, y=.data$ReadID)) + 
      geom_tile(aes(fill=.data$Methylation), height=1, width=1) +
      facet_wrap(~.data$Sample, scales = "free_y", dir = 'v') +
      ylab("") +
      xlab("") +
      scale_discrete_manual(aesthetics = "fill", values = c("black", "grey")) +
      theme_classic() +
      theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.line = element_blank()) -> pl
  }
  
  return(pl)
  
}

#' Single TF state quantification bar
#'
#' @param SortedReads Sorted reads as returned by SortReadsBySingleTF
#'
SingleTFStateQuantificationPlot = function(SortedReads){
  
  StateQuantificationPlot(SortedReads = SortedReads, states = SingleTFStates())
  
}

#' TF pair state quantification bar
#'
#' @param SortedReads Sorted reads as returned by SortReadsByTFCluster
#'
TFPairStateQuantificationPlot = function(SortedReads){
  
  StateQuantificationPlot(SortedReads = SortedReads, states = TFPairStates())
  
}

#' Plot SMF data at single site
#'
#' @param Methylation Context methylation object as returned by CallContextMethylation function
#' @param RegionOfInterest GRanges interval to plot
#' @param ShowContext TRUE or FALSE (default). Causes the genomic context of the plotted cytosines to be displayed as the dot shape
#' @param TFBSs GRanges object of transcription factor binding sites to include in the plot. Assumed to be already subset. Also assumed that the tf names are under the column "TF"
#' @param SNPs GRanges object of SNPs to visualize. Assumed to be already subset. Assumed to have the reference and alternative sequences respectively under the columns "R" and "A"
#' @param SortingBins GRanges object of sorting bins (absolute) coordinate to visualize
#' @param SortedReads Defaults to NULL, in which case will plot unsorted reads. Sorted reads object as returned by SortReads function or "HC" to perform hierarchical clustering
#' @param sorting.strategy One of "classical" (default), "custom", "hierarchical.clustering" or "None". Determines how to display reads. For details check documentation from PlotSM function.
#'
#' @importFrom grDevices dev.list dev.off pdf
#' @importFrom patchwork plot_layout
#'
#' @export
#'
#' @examples
#' 
#' library(GenomicRanges)
#'
#' RegionOfInterest = GRanges("chr12", IRanges(20464551, 20465050))
#' Methylation = qs::qread(system.file("extdata", "Methylation_3.qs", 
#' package="SingleMoleculeFootprinting"))
#' TFBSs = qs::qread(system.file("extdata", "TFBSs_3.qs", package="SingleMoleculeFootprinting"))
#' SortedReads = SortReadsBySingleTF(MethSM = Methylation[[2]], TFBS = TFBSs)
#'
#' PlotSingleSiteSMF(Methylation = Methylation,
#'                   RegionOfInterest = RegionOfInterest,
#'                   SortedReads = SortedReads,
#'                   TFBSs = TFBSs)
#'
PlotSingleSiteSMF = function(Methylation, RegionOfInterest, ShowContext=FALSE, TFBSs=NULL, SNPs=NULL, SortingBins=NULL, SortedReads=NULL, sorting.strategy = "None"){

  message("Producing average SMF plot")
  PlotAvgSMF(MethGR = Methylation[[1]],
             MethSM = Methylation[[2]],
             RegionOfInterest = RegionOfInterest,
             SortedReads = SortedReads,
             ShowContext = ShowContext,
             TFBSs = TFBSs,
             SNPs = SNPs,
             SortingBins = SortingBins) -> Avg_pl

  message("Producing Single Molecule stacks")
  PlotSM(MethSM = Methylation[[2]], RegionOfInterest = RegionOfInterest, SortedReads = SortedReads, sorting.strategy = sorting.strategy) -> SM_pl

  # State quantification plot
  if(is.list(SortedReads)){
    message("Producing state quantification plots")
    
    if(all(nchar(names(SortedReads[[1]])) == 3)){
      states = SingleTFStates()
      StateQuantificationPlot(SortedReads = SortedReads, states = states) -> StateQuant_pl
    } else if(all(nchar(names(SortedReads[[1]])) == 4)){
      states = TFPairStates()
      StateQuantificationPlot(SortedReads = SortedReads, states = states) -> StateQuant_pl
    } else {
      message("No states recognised, skipping StateQuantificationPlot")
      StateQuant_pl = NULL
    }
    
  } else {
    StateQuant_pl = NULL
  }
  
  message("Combining plots")
  layout <- "
  #A
  CB
  "
  Avg_pl + SM_pl + StateQuant_pl +
    patchwork::plot_layout(ncol = 2, design = layout, widths = c(0.25, 1), heights = c(1, 0.8), guides = "collect") -> FinalPlot
  
  return(FinalPlot)

}

