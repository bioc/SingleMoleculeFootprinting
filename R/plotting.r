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
#' Qinput = system.file("extdata", "QuasR_input_pairs.txt", package = "SingleMoleculeFootprinting", mustWork = T)
#' MySample = suppressMessages(readr::read_delim(Qinput, delim = "\t")[[2]])
#' Region_of_interest = GRanges(seqnames = "chr6", ranges = IRanges(start = 88106000, end = 88106500), strand = "*")
#' Methylation = CallContextMethylation(sampleSheet = Qinput,
#'                                      sample = MySample,
#'                                      genome = BSgenome.Mmusculus.UCSC.mm10,
#'                                      range = Region_of_interest,
#'                                      coverage = 20,
#'                                      ConvRate.thr = 0.2)
#' TFBSs = GenomicRanges::GRanges("chr6", IRanges(c(88106253), c(88106263)), strand = "-")
#' elementMetadata(TFBSs)$name = c("NRF1")
#' names(TFBSs) = c(paste0("TFBS_", c(4305216)))
#'
#' PlotAvgSMF(MethGR = Methylation[[1]], Region_of_interest = Region_of_interest, TFBSs = TFBSs)
#'
PlotAvgSMF = function(MethGR, MethSM=NULL, RegionOfInterest, SortedReads=NULL, ShowContext=FALSE, TFBSs=NULL, SNPs=NULL, SortingBins=NULL){

  # Prepare SMF data
  MethGR %>%
    as_tibble() %>%
    select(-grep("_Coverage$", colnames(.)), -end, -width, -strand) %>%
    gather(sample, MethRate, -seqnames, -start, -GenomicContext) %>%
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
      select(start, end, TF) -> TFBS_PlottingDF
  }

  # Prepare SNPs
  if(!is.null(SNPs)){
    SNPs %>%
      as_tibble() %>%
      rowwise() %>%
      mutate(width = ifelse(max(nchar(R), nchar(A)) == 1, 3, max(nchar(R), nchar(A)))) %>%
      select(start, width) %>%
      mutate(y_coord = -0.13) -> SNPs_PlottingDF
  }

  # Prepare SortingBins
  if(!is.null(SortingBins)){
    SortingBins %>%
      as.data.frame() %>%
      select(start, end) -> Bins_PlottingDF
  }

  PlottingDF %>%
    ggplot(aes(x=start, y=1-MethRate, color=sample)) +
    geom_line() +
    {if(ShowContext){geom_point(aes(shape=GenomicContext))}else{geom_point()}} +
    {if(!is.null(TFBSs)){geom_rect(TFBS_PlottingDF, mapping = aes(xmin=start, xmax=end, ymin=-0.09, ymax=-0.04), inherit.aes = FALSE)}} +
    {if(!is.null(TFBSs)){ggrepel::geom_text_repel(TFBS_PlottingDF, mapping = aes(x=start+((end-start)/2), y=-0.02, label=TF), min.segment.length = .1, max.overlaps = 1e+05, inherit.aes = FALSE)}} + #, size=4.5
    {if(!is.null(SNPs)){geom_tile(SNPs_PlottingDF, mapping = aes(x=start, y=y_coord, width=width), color = ColorsToUse[2], fill = ColorsToUse[2], height = 0.05, inherit.aes = FALSE)}} +
    {if(!is.null(SortingBins)){geom_rect(Bins_PlottingDF, mapping = aes(xmin=start, xmax=end, ymin=-0.02, ymax=0), color="black", fill="white", inherit.aes = FALSE)}} +
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
#'
#' @export
#'
#' @examples
#'
#' Qinput = system.file("extdata", "QuasR_input_pairs.txt", package = "SingleMoleculeFootprinting", mustWork = T)
#' MySample = suppressMessages(readr::read_delim(Qinput, delim = "\t")[[2]])
#' Region_of_interest = GRanges(seqnames = "chr6", ranges = IRanges(start = 88106000, end = 88106500), strand = "*")
#' Methylation = CallContextMethylation(sampleSheet = Qinput,
#'                                      sample = MySample,
#'                                      genome = BSgenome.Mmusculus.UCSC.mm10,
#'                                      range = Region_of_interest,
#'                                      coverage = 20,
#'                                      ConvRate.thr = 0.2)
#'
#' PlotSingleMoleculeStack(MethSM = Methylation[[2]], RegionOfInterest = Region_of_interest)
#'
PlotSingleMoleculeStack = function(MethSM, RegionOfInterest){
  
  Reduce(rbind, lapply(seq_along(MethSM), function(i){
    MethSM[[i]] %>%
      as.matrix() %>%
      as.data.frame() %>%
      rownames_to_column(var = "ReadID") %>%
      mutate(Sample = names(MethSM)[i]) %>%
      gather(Coordinate, Methylation, -ReadID, -Sample) %>%
      mutate(Methylation = ifelse(Methylation == 0, NA, Methylation-1))
  })) %>%
    na.omit() %>%
    mutate(Methylation = as.factor(Methylation), Coordinate = as.numeric(Coordinate)) -> PlottingDF
  PlottingDF$ReadID = factor(PlottingDF$ReadID, levels = Reduce(c, lapply(MethSM, rownames)))
  
  PlottingDF %>%
    group_by(Sample) %>%
    summarise(NrReads = length(unique(ReadID))) %>%
    ungroup() %>%
    mutate(Label = paste0(Sample, " (", NrReads, " reads)")) %>%
    select(Sample, Label) -> Labels
  names(Labels$Label) = Labels$Sample

  PlottingDF %>%
    ggplot(aes(x=Coordinate, y=ReadID)) + 
    geom_tile(aes(fill=Methylation), height=1, width=5) +
    facet_wrap(~Sample, scales = "free_y", dir = 'v', 
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
#' Qinput = system.file("extdata", "QuasR_input_pairs.txt", package = "SingleMoleculeFootprinting", mustWork = T)
#' MySample = suppressMessages(readr::read_delim(Qinput, delim = "\t")[[2]])
#' Region_of_interest = GRanges(seqnames = "chr6", ranges = IRanges(start = 88106000, end = 88106500), strand = "*")
#' Methylation = CallContextMethylation(sampleSheet = Qinput,
#'                                      sample = MySample,
#'                                      genome = BSgenome.Mmusculus.UCSC.mm10,
#'                                      RegionOfInterest = Region_of_interest,
#'                                      coverage = 20,
#'                                      ConvRate.thr = 0.2)
#'
#'  PlotSM(MethSM = Methylation[[2]], RegionOfInterest = Region_of_interest)
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
#'
#' @export
#'
#' @examples
#'
#' Qinput = system.file("extdata", "QuasR_input_pairs.txt", package = "SingleMoleculeFootprinting", mustWork = T)
#' MySample = suppressMessages(readr::read_delim(Qinput, delim = "\t")[[2]])
#' Region_of_interest = GRanges(seqnames = "chr6", ranges = IRanges(start = 88106000, end = 88106500), strand = "*")
#' Methylation = CallContextMethylation(sampleSheet = Qinput,
#'                                      sample = MySample,
#'                                      genome = BSgenome.Mmusculus.UCSC.mm10,
#'                                      range = Region_of_interest,
#'                                      coverage = 20,
#'                                      ConvRate.thr = 0.2)
#' TFBSs = GenomicRanges::GRanges("chr6", IRanges(c(88106253), c(88106263)), strand = "-")
#' elementMetadata(TFBSs)$name = c("NRF1")
#' names(TFBSs) = c(paste0("TFBS_", c(4305216)))
#' SortedReads = SortReadsByTFCluster(MethSM = Methylation[[2]], TFBSs = TFBSs)
#'
#' StateQuantificationPlot(SortedReads = SortedReads, states = TFPairStates())
#'
StateQuantificationPlot = function(SortedReads, states){
  
  if(states[[1]] == "101"){
    states_df = unique(gather(as_tibble(states), "State", "Pattern")) %>%
      mutate(State = factor(State, levels = c("bound", "accessible", "closed", "unassigned")))
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
    arrange(desc(State)) %>%
    mutate(ReadID = factor(ReadID, levels = ReadID)) -> PlottingDF
  
  if(states[[1]] == "101"){
    PlottingDF %>%
      ggplot(aes(x=1, y=ReadID, fill=State)) + 
      geom_tile() +
      facet_wrap(~Sample, scales = "free_y", dir = 'v') +
      ylab("") +
      xlab("") +
      scale_fill_manual(breaks = names(states), values = c("#984EA3", "#4DAF4A", "#377EB8", "#999999")) + 
      theme_classic() +
      theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.line = element_blank()) -> pl
  } else if (states[[1]] == "1001") {
    PlottingDF %>%
      separate(Pattern, into = c(paste0("Bin", seq(unique(nchar(unlist(states)))))), sep = "(?<=.)", extra = 'drop') %>%
      gather(Bin, Methylation, -ReadID, -Sample, -State) %>%
      ggplot(aes(x=Bin, y=ReadID)) + 
      geom_tile(aes(fill=Methylation), height=1, width=1) +
      facet_wrap(~Sample, scales = "free_y", dir = 'v') +
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
#' Qinput = system.file("extdata", "QuasR_input_pairs.txt", package = "SingleMoleculeFootprinting", mustWork = T)
#' MySample = suppressMessages(readr::read_delim(Qinput, delim = "\t")[[2]])
#' Region_of_interest = GRanges(seqnames = "chr6", ranges = IRanges(start = 88106000, end = 88106500), strand = "*")
#' Methylation = CallContextMethylation(sampleSheet = Qinput,
#'                                      sample = MySample,
#'                                      genome = BSgenome.Mmusculus.UCSC.mm10,
#'                                      range = Region_of_interest,
#'                                      coverage = 20,
#'                                      ConvRate.thr = 0.2)
#' TFBSs = GenomicRanges::GRanges("chr6", IRanges(c(88106253), c(88106263)), strand = "-")
#' elementMetadata(TFBSs)$name = c("NRF1")
#' names(TFBSs) = c(paste0("TFBS_", c(4305216)))
#' SortedReads = SortReadsByTFCluster(MethSM = Methylation[[2]], TFBSs = TFBSs)
#'
#' PlotSingleSiteSMF(ContextMethylation = Methylation,
#'                   sample = MySample,
#'                   range = Region_of_interest,
#'                   SortedReads = SortedReads,
#'                   TFBSs = TFBSs,
#'                   saveAs = NULL)
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
    StateQuantificationPlot(SortedReads = SortedReads) -> StateQuant_pl
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

