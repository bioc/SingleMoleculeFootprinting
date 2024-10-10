# A-posteriori annotate TF footprints by 
#     (1) Asking whether we have an already annotated TFBS from the list passed upstream by the user
#     (2) Scoring the fooprinted sequence for all known PWMs and picking the highest-scoring TF
#     (2.note) No in vivo occupancy data is used to verify the motif prediction since doing so
#              does not return predictions (even at lenient parameters)

#' This just crosses info between the automatically discovered TF footprints and the user-provided TFBSs on which the analysis has been centered
#' 
#' @param footprints.df as coming frmo detect.footprint.wrapper
#' @param chromosome
#' @param TFBSs.overlapping.locus already annotated TFBSs passed upstream by the user
#' @param assign.TF.identity whether to actually assign the TF id to footprints or to just use the info from TFBSs.overlapping.locus to manipulate footprints.df. Defaults to TRUE
#' 
assign.passed.TFBSs.to.footprints = function(footprints.df, chromosome, TFBSs.overlapping.locus, assign.TF.identity=TRUE){
  
  footprints.df %>%
    arrange(start, desc(width)) %>%
    mutate(seqnames = chromosome) %>%
    GRanges() -> footprints.gr
  
  footprints.gr$TF = NA
  footprints.gr$TF.name = NA
  if(assign.TF.identity){
    Overlaps = findOverlaps(footprints.gr, IRanges::resize(TFBSs.overlapping.locus, 1, "center")) # N.b.: by resizing to 1, we require that at least half the TFBSs falls into the footprint
    footprints.gr$TF[unique(queryHits(Overlaps))] = split(as.character(TFBSs.overlapping.locus$TF)[subjectHits(Overlaps)], queryHits(Overlaps))
    footprints.gr$TF.name[unique(queryHits(Overlaps))] = split(as.character(TFBSs.overlapping.locus$absolute.idx)[subjectHits(Overlaps)], queryHits(Overlaps))
    footprints.gr[footprints.gr$biological.state != "TF"]$TF = NA # remove TF from non-TF footprints
    footprints.gr[footprints.gr$biological.state != "TF"]$TF.name = NA
  }
  
  footprints.df = as.data.frame(footprints.gr)
  
  return(footprints.df)
  
}

#' Map user-provided PWMs to so-far-orphan footprints
#' 
#' N.b. assuming the user provided upstream a list of TFBSs with in-vivo evidence of occupancy (e.g. ChIPseq, ChIPexo etc..)
#' the TFBSs obtained with this function are to be considered less trustworthy given the lack of in-vivo evidence
#'
#' @param footprints.df as coming from detect.footprint.wrapper (or assign.passed.TFBSs.to.footprints)
#' @param PWM.list List of PWMs. If set to NULL, JASPAR2020 PWMs will be used, see prepare.JASPAR.PWMs for details.
#' @param genome
#' @param cores number of cores to use to scan PWMs, defaults to 1.
#'
#' @import TFBSTools
#' @import JASPAR2020
#' @import dplyr
#' @importFrom plyranges reduce_ranges
#' 
library(TFBSTools)
annotate.orphan.footprints = function(footprints.df, params){
  
  # Unpack parameters ####
  PWM.list = params$annotate.orphan.footprints[["PWM.list"]]
  genome = params$annotate.orphan.footprints[["genome"]]
  cores = params$annotate.orphan.footprints[["cores"]]
  ########################
  
  if(is.null(footprints.df$TF)){
    footprints.df$TF = NA # if empty, initialize
  }
  
  footprints.df[is.na(footprints.df$TF) & footprints.df$biological.state == "TF",] %>%
    select(seqnames, start, end, width, strand) %>%
    distinct() %>%
    GRanges() -> orphan.footprints.gr
  message(paste0("Attempting to match ", length(orphan.footprints.gr), " orphan footprints with a TF"))
  
  if (is.null(PWM.list)){
    stop("Please pass PWM.list")
  }
  
  seq2query = Biostrings::getSeq(genome, orphan.footprints.gr)
  sitesetList = searchSeq(PWM.list, seq2query, seqname=seqnames(orphan.footprints.gr), min.score=0.85, strand="*", mc.cores = cores)
  Reduce(rbind, lapply(sitesetList, function(x){TFBSTools::as.data.frame(x)})) -> match.df
  message(paste0(nrow(match.df), " matches found"))
  
  match.df %>%
    arrange(seqnames, desc(absScore)) %>%
    group_by(seqnames) %>%
    slice_max(order_by = absScore, n = 3) %>%
    ungroup() %>%
    select(seqnames, absScore, TF) %>%
    mutate(seqnames = as.integer(seqnames))-> top.scores
  
  # gather data for PWM scanning diagnostic plot
  orphan.footprints.gr$orphan.fp.idx = as.integer(seq_along(orphan.footprints.gr))
  orphan.footprints.gr$TF = NULL
  PWM.scanning.results = dplyr::full_join(x = as_tibble(orphan.footprints.gr), y = top.scores, by = c("orphan.fp.idx" = "seqnames"))
  
  # annotate orphan footprints 
  missing.TFs = rep(
    slice_max(group_by(PWM.scanning.results, orphan.fp.idx), order_by = "asbScore", n = 1)$TF,
    table(paste(footprints.df[is.na(footprints.df$TF) & footprints.df$biological.state == "TF",]$start, 
                footprints.df[is.na(footprints.df$TF) & footprints.df$biological.state == "TF",]$end, "-"))
  )
  footprints.df$TF.annotation = NA
  footprints.df$TF.annotation[!is.na(footprints.df$TF)] = "user.provided"
  footprints.df[is.na(footprints.df$TF) & footprints.df$biological.state == "TF",]$TF = missing.TFs
  footprints.df$TF.annotation[footprints.df$biological.state == "TF" & is.na(footprints.df$TF.annotation)] = "rescanned"
  
  return(list(footprints.df = footprints.df, PWM.scanning.results = PWM.scanning.results))
  
}

#' N.b. atm will annotate with archetypes also the user.passed TF annotations
#' 
#' @param footprints.df has to come from annotate.orphan.footprints
#' @param archetypes
#' 
reduce.to.archetypes = function(footprints.df, archetypes){
  
  query.TFs = footprints.df$TF[!is.na(footprints.df$TF)]
  
  all.components = strsplit(archetypes$cluster.components, ",")
  names(all.components) = paste0(archetypes$central.motif, "_")
  all.components.split = unlist(all.components, use.names = TRUE)
  
  footprints.df$archetype.annotation = NA
  footprints.df[!is.na(footprints.df$TF),][query.TFs %in% all.components.split,]$archetype.annotation = 
    gsub("_.*$", "", names(unlist(lapply(query.TFs[query.TFs %in% all.components.split], function(tf){which(all.components.split == tf)[1]}))))
  
  return(footprints.df)
  
}

#' Overall wrapper for footprint annotation process
#' 
#' @param footprints.df as coming from detect.footprint.wrapper
#' @param chromosome chromosome of current Region of interest. Needed for annotation
#' @param TFBSs.overlapping.locus already annotated TFBSs passed upstream by the user
#' @param RNAseq RNAseq-derived count table (e.g. DESeq2) to interpolate and ask whether TF identifyied is expressed (hence corroborating annotation)
#' @param params
#' 
Annotate.footprints.wrapper = function(footprints.df, chromosome, TFBSs.overlapping.locus, RNAseq=NULL, params){
  
  # Unpack parameters ####
  scan.orphan.footprints = params$annotate.orphan.footprints[["scan.orphan.footprints"]]
  use.passed.TFBSs.for.annotation = params$annotate.orphan.footprints[["use.passed.TFBSs.for.annotation"]]
  archetypes = params$annotate.orphan.footprints[["archetypes"]]
  ########################
  
  message("### Annotating footprints ###")
  if(use.passed.TFBSs.for.annotation){
    message("# 1. Passing user-provided TFBSs to annotate footprints")
    footprints.df = assign.passed.TFBSs.to.footprints(footprints.df = footprints.df, chromosome = chromosome, TFBSs.overlapping.locus = TFBSs.overlapping.locus, assign.TF.identity = TRUE)
    TFBSs.overlapping.locus$TF.annotation = rep("user.provided", length(TFBSs.overlapping.locus))
  } else {
    footprints.df = assign.passed.TFBSs.to.footprints(footprints.df = footprints.df, chromosome = chromosome, TFBSs.overlapping.locus = TFBSs.overlapping.locus, assign.TF.identity = FALSE)
    TFBSs.overlapping.locus = TFBSs.overlapping.locus[0]
  }
  
  if(sum(is.na(footprints.df[footprints.df$biological.state == "TF",]$TF)) == 0){
    message("# 2. All TF footprints assigned...skipping next step")
  } else if (scan.orphan.footprints == FALSE){
    message("# 2. skipping rescanning of orphan footprints")
  } else {
    message("# 2. Scanning DNA sequence of thus-far orphan footprints for PWMs")
    annotation.results = annotate.orphan.footprints(footprints.df = footprints.df, params=params)
    footprints.df = annotation.results$footprints.df
    PWM.scanning.results = annotation.results$PWM.scanning.results
    TFBSs.overlapping.locus = c(
      TFBSs.overlapping.locus,
      PWM.scanning.results %>% 
        group_by(orphan.fp.idx) %>%
        slice_max(order_by = absScore, n = 1) %>%
        ungroup() %>%
        mutate(TF.annotation = "rescanned") %>%
        select(seqnames, start, end, width, strand, TF, TF.annotation) %>%
        GRanges()
    )
    
    if(!is.null(archetypes)){
      message("# 2.1. Reducing annotated TFs to archetypes")
      footprints.df = reduce.to.archetypes(footprints.df = footprints.df, archetypes = archetypes)
    }
    
  }
  
  if(!is.null(RNAseq)){
    warning("IMPLEMENT ME PLEASE")
  }
  
  message("# 3. Gathering results")
  if(!is.null(TFBSs.overlapping.locus)){
    TFBSs.overlapping.locus = TFBSs.overlapping.locus %>% plyranges::arrange(start, end)
  }
  
  results = list(footprints.df = footprints.df, TFBSs.overlapping.locus = TFBSs.overlapping.locus)
  
  return(results)
  
}


