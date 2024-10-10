Load.data = function(data.type = NULL, pool.replicates = FALSE){
  
  current.variables = ls(name = sys.frame(0))
  
  if(isEmpty(grep("^sampleSheet$", current.variables))){
    if(data.type == "F1_bait.capture"){
      if(pool.replicates){sampleSheet <<- "/g/krebs/barzaghi/HTS/SMF/MM/QuasR_input_all_bait_capture_data_pooledReplicates.txt"} 
      else if (!pool.replicates){sampleSheet <<- "/g/krebs/barzaghi/HTS/SMF/MM/QuasR_input_all_bait_capture_data.txt"}
    } else if (data.type == "F1_amplicon"){
      if(pool.replicates){sampleSheet <<- "/g/krebs/barzaghi/HTS/SMF/MM/QuasR_input_all_amplicon_data.txt"} 
      else if (!pool.replicates){stop("We didn't make this one")}
    } else if(data.type == "WT_bait.capture"){
      if(pool.replicates){sampleSheet <<- "/g/krebs/barzaghi/HTS/SMF/MM/QuasR_input_files/QuasR_input_AllCanWGpooled_dprm_DE_only.txt"} 
      else if (!pool.replicates){stop("We didn't make this one")}
    } else if (data.type == "WT_amplicon"){
      if(pool.replicates){sampleSheet <<- "/g/krebs/barzaghi/HTS/SMF/MM/QuasR_input_files/Can_amplicons_NRF1KD_QuasR_input.txt"}
      else if (!pool.replicates){stop("We didn't make this one")}
    } else if (data.type == "Rest_ko_amplicon"){
      if(pool.replicates){sampleSheet <<- "/g/krebs/barzaghi/HTS/SMF/MM/AMP_REST_NO/samplesheet.txt"}
      else if (!pool.replicates){stop("We didn't make this one")}
    } else if (data.type == "Nrf1_kd_amplicon"){
      if(pool.replicates){sampleSheet <<- "/g/krebs/barzaghi/HTS/SMF/MM/QuasR_input_files/Can_amplicons_NRF1KD_QuasR_input.txt"}
      else if (!pool.replicates){stop("We didn't make this one")}
    } else if (data.type == "Sox2_kd_bait.capture"){
      if(pool.replicates){sampleSheet <<- "/g/krebs/barzaghi/HTS/SMF/MM/2024-05-22-2222T2JNX/Qinput_Sox2_dTAG.txt"}
      else if (!pool.replicates){stop("We didn't make this one")}
    } else if (data.type == "Oct4_kd_bait.capture"){
      if(pool.replicates){sampleSheet <<- "/g/krebs/barzaghi/HTS/SMF/MM/2024-05-22-2222T2JNX/Qinput_Oct4_dTAG.txt"}
      else if (!pool.replicates){stop("We didn't make this one")}
    } else if (data.type == "remodeller_bait.capture"){
      if(pool.replicates){sampleSheet <<- "/g/krebs/barzaghi/HTS/SMF/MM/2024-05-22-2222T2JNX/Qinput_remodeller.txt"}
      else if (!pool.replicates){stop("We didn't make this one")}
    } else if (data.type == "PRA_amplicon"){
      if(pool.replicates){sampleSheet <<- "/g/krebs/barzaghi/analyses/05.04.23_Unsupervised_on_PRA/Qinput_all_insertion_data_RepsPooled.txt"}
      else if (!pool.replicates){sampleSheet <<- "/g/krebs/barzaghi/analyses/05.04.23_Unsupervised_on_PRA/Qinput_all_insertion_data.txt"}
    } else if (data.type == "PRA_mutant_amplicon"){
      if(pool.replicates){sampleSheet <<- "/g/krebs/barzaghi/analyses/08.02.24_RMCE_mutants/samplesheet_merged.txt"}
      else if (!pool.replicates){sampleSheet <<- "/g/krebs/barzaghi/analyses/08.02.24_RMCE_mutants/samplesheet.txt"}
    }
  }
    
  if(isEmpty(grep("^MySamples$", current.variables))){
    if(data.type %in% c("F1_bait.capture", "Sox2_kd_bait.capture", "Oct4_kd_bait.capture")){
      MySamples <<- list(NO = grep("NO", readr::read_delim(sampleSheet, "\t")$SampleName, value = TRUE) %>% unique(), DE = grep("DE", readr::read_delim(sampleSheet, "\t")$SampleName, value = TRUE) %>% unique())
    } else if (data.type %in% c("F1_amplicon", "WT_bait.capture", "WT_amplicon", "Nrf1_kd_amplicon", "remodeller_bait.capture", "PRA_amplicon", "PRA_mutant_amplicon")){
      MySamples <<- list(NO = NULL, DE = readr::read_delim(sampleSheet, delim = "\t", col_names = TRUE)$SampleName %>% unique())
    } else if (data.type == "Rest_ko_amplicon"){
      MySamples <<- list(NO = readr::read_delim(sampleSheet, delim = "\t", col_names = TRUE)$SampleName %>% unique(), DE = NULL)
    }
  }
  if(isEmpty(grep("^TFBSs$", current.variables))){
    if(data.type == "PRA_mutant_amplicon"){
      message("loading TFBSs")
      TFBSs <<- qs::qread("/g/krebs/barzaghi/analyses/08.02.24_RMCE_mutants/TFBSs.qs")
    } else {
      message("loading TFBSs")
      TFBSs <<- Load.TFBSs(kind = "Barzaghi", GenVar.info = "GenVarChange")
    }}
  if (data.type %in% c("F1_bait.capture", "F1_amplicon")){
      if(isEmpty(grep("^SNPs$", current.variables))){
      message("loading SNPs")
      SNPs <<- LoadSNPs()}
    if(isEmpty(grep("^CytosinesToMask$", current.variables))){
      message("loading CytosinesToMask")
      CytosinesToMask <<- Load.CytosinesToMask()}
  } else if (data.type %in% c("WT_bait.capture", "WT_amplicon", "Rest_ko_amplicon", "Nrf1_kd_amplicon", "Sox2_kd_bait.capture", "Oct4_kd_bait.capture", "remodeller_bait.capture", "PRA_amplicon", "PRA_mutant_amplicon")){
      SNPs <<- NULL
      CytosinesToMask <<- NULL
    }
  if(isEmpty(grep("^GenomicTiles$", current.variables))){
    message("Loading genomic tiles")
    if(data.type %in% c("F1_bait.capture", "WT_bait.capture", "Sox2_kd_bait.capture", "Oct4_kd_bait.capture", "remodeller_bait.capture")){
      GenomicTiles <<- Load.GenomicTiles(tiles.width = 80, GeneAssociation = FALSE, CGI = TRUE)
    } else if (data.type == "F1_amplicon"){
      GenomicTiles <<- c(
        readRDS("/g/krebs/barzaghi/analyses/F1_Amplicons/Run1/final.amplicon.ranges.rds"),
        readRDS("/g/krebs/barzaghi/analyses/F1_Amplicons/Run2/final.amplicon.ranges.rds")) %>% sort()
      names(GenomicTiles) = paste0("AMP_", seq_along(GenomicTiles))
      GenomicTiles <<- GenomicTiles
    } else if (data.type %in% c("WT_amplicon", "Rest_ko_amplicon", "Nrf1_kd_amplicon")){
      GenomicTiles <<- Load.Sonmezer.amplicon.GRanges()
    } else if (data.type == "PRA_amplicon"){
      GenomicTiles <<- qs::qread("/g/krebs/barzaghi/analyses/05.04.23_Unsupervised_on_PRA/CRE_autonomy_final_ranges.qs")
    } else if (data.type == "PRA_mutant_amplicon"){
      GenomicTiles <<- qs::qread("/g/krebs/barzaghi/analyses/08.02.24_RMCE_mutants/MutantsBatch1VB_ranges.qs")
    }
  }
  
}

plotting.color.assignment = function(current.sample){
  if(current.sample %in% c("SMF_MM_TKO_DE_", "Rest_ko", "amplicon_DE_data", "SMF_MM_TKO_DE_ectopic", "MUT_DE_") | grepl("MUTR[1-3][a-c]", current.sample)){
    plotting.colors = "black"
  } else if(current.sample %in% c("NRF1_KD_DE_", "Sox2_kd_NO_", "Oct4_kd_NO_", "remodeller")){
    plotting.colors = c("black", "red")
  } else if(current.sample == "STKO"){
    plotting.colors = c("black", "grey45")
  } else if(current.sample == "CTKO"){
    plotting.colors = c("black", "sienna")
  } else {
    stop("current.sample not recognised...plotting.color assignment failed")
  }
  return(plotting.colors)
}

#'
polish.bulk.smf.plot = function(footprint.charter.results, data.type, current.sample, RegionOfInterest, plotting.TFBSs){
  
  footprint.charter.results$SingleSite.plot$data$.GenomicContext = footprint.charter.results$SingleSite.plot$data$GenomicContext
  footprint.charter.results$SingleSite.plot$data$GenomicContext = ""
  
  plotting.TFBSs.layer.data = plotting.TFBSs %>% data.frame() %>% dplyr::select(start, end, TF)
  footprint.charter.results$SingleSite.plot$layers[[4]]$data = plotting.TFBSs.layer.data
  footprint.charter.results$SingleSite.plot$layers[[3]]$data = plotting.TFBSs.layer.data
  
  x.axis.breaks = c(start(resize(footprint.charter.results$RegionOfInterest, width = 500, fix = "center")), end(resize(footprint.charter.results$RegionOfInterest, width = 500, fix = "center")))
  
  footprint.charter.results$SingleSite.plot +
    geom_point(size = 2.5) +
    geom_line(linewidth = 1.25) +
    scale_color_manual(values = plotting.color.assignment(current.sample)) +
    scale_fill_manual(values = plotting.color.assignment(current.sample)) +
    scale_x_continuous(limits = x.axis.breaks, breaks = x.axis.breaks, labels = format(x.axis.breaks, nsmall=1, big.mark=",")) +
    scale_y_continuous(limits = c(-ifelse(data.type %in% c("F1_bait.capture", "F1_amplicon"), .17, .1),1), breaks = c(0,1), labels = function(x) format(100*x,digits = 2)) +
    guides(shape = "none") +
    ggtitle(NULL) +
    xlab(as.character(seqnames(RegionOfInterest))) +
    ylab("SMF (1 - meth %)") +
    theme_bw() +
    theme(text = element_text(size = 18), axis.title.y = element_text(vjust = -5), axis.title.x = element_text(vjust = 6), 
          panel.border = element_rect(linewidth = 1, fill = "transparent"))
  
}

#'
polish.sm.stacks = function(footprint.charter.results, current.sample, partition.collapsing.dict){
  
  lapply(seq_along(current.sample), function(i){
    if(!is.null(partition.collapsing.dict)){
      plotting.reads = footprint.charter.results$pam.res[names(footprint.charter.results$pam.res) %in% footprint.charter.results$SM.plot[[current.sample[i]]]$data$ReadID]
      plotting.reads = lapply(partition.collapsing.dict, function(x){names(plotting.reads)[plotting.reads %in% x]})
      left_join(
        footprint.charter.results$SM.plot[current.sample][[i]]$data,
        rownames_to_column(dplyr::rename(as.data.frame(footprint.charter.results$pam.res), "cluster.nr" = `footprint.charter.results$pam.res`), var = "ReadID"), 
        by = "ReadID"
      ) %>% mutate(ReadID = factor(ReadID, levels = rev(as.character(unlist(plotting.reads))))) -> PlottingDF
    } else {
      left_join(
        footprint.charter.results$SM.plot[current.sample][[i]]$data,
        rownames_to_column(dplyr::rename(as.data.frame(footprint.charter.results$pam.res), "cluster.nr" = `footprint.charter.results$pam.res`), var = "ReadID"), 
        by = "ReadID"
      ) -> PlottingDF
    }
    nr.molecules = length(unique(PlottingDF$ReadID))
    x.axis.breaks = c(start(resize(footprint.charter.results$RegionOfInterest, width = 500, fix = "center")), end(resize(footprint.charter.results$RegionOfInterest, width = 500, fix = "center")))
    PlottingDF %>%
      ggplot(aes(x = Coordinate, y = ReadID)) + 
      geom_tile(aes(fill = Methylation), height = 1, width = 5) + 
      {if(is.null(partition.collapsing.dict)){facet_grid(factor(cluster.nr, levels = suppressWarnings(arrange.clusters.by.footprints(footprint.charter.results$footprints.df)))~Sample, scales = "free", space = "free_y")}} +
      ylab(paste0(nr.molecules, " molecules")) +
      scale_x_continuous(limits = x.axis.breaks, breaks = x.axis.breaks, labels = format(x.axis.breaks, nsmall=1, big.mark=",")) +
      scale_discrete_manual(aesthetics = "fill", values = c("black", "grey")) + 
      theme_classic() + 
      xlab(as.character(seqnames(footprint.charter.results$RegionOfInterest))) + 
      theme(text = element_text(size = 18), axis.title.y = element_text(vjust = 0), axis.title.x = element_text(vjust = 5),
            axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.border = element_rect(linewidth = 1, fill = "transparent"))
  }) -> SM.plots
  names(SM.plots) = gsub("_.*Barzaghi|_MethRate", "", current.sample)
  
  return(SM.plots)
  
}

#'
klee.plot = function(footprint.charter.results, current.sample, partition.collapsing.dict){
  
  x.axis.breaks = c(start(resize(footprint.charter.results$RegionOfInterest, width = 500, fix = "center")), end(resize(footprint.charter.results$RegionOfInterest, width = 500, fix = "center")))
  
  footprint.charter.results$footprints.df %>%
    rowwise() %>%
    mutate(start = list(seq(start, end))) %>%
    ungroup() %>%
    dplyr::select(start, biological.state, cluster.nr, cluster.coverage) -> fp.df
  
  lapply(current.sample, function(current.sample.internal){
    
    left_join(
      fp.df %>% dplyr::select(-cluster.coverage),
      footprint.charter.results$cluster.coverage.df %>%
        filter(Sample == current.sample.internal) %>%
        dplyr::select(-TFBS.cluster),
      by = c("cluster.nr")
    ) -> fp.df
    
    if(!is.null(partition.collapsing.dict)){
      fp.df %<>% 
        mutate(cluster.nr = factor(cluster.nr, levels = rev(unlist(partition.collapsing.dict, use.names = FALSE)))) %>%
        arrange(cluster.nr)
    } else {
      fp.df %<>% 
        mutate(cluster.nr = factor(cluster.nr, levels = suppressWarnings(arrange.clusters.by.footprints(footprint.charter.results$footprints.df)))) %>%
        arrange(cluster.nr)
    }
    
    fp.df %>%
      dplyr::select(cluster.nr, read.count, Sample) %>%
      distinct() %>%
      group_by(Sample) %>%
      mutate(read.count.cumsum = cumsum(read.count), d = read.count.cumsum - read.count + 1) %>%
      rowwise() %>%
      mutate(read.idx = list(seq(d,read.count.cumsum))) %>%
      ungroup() %>%
      dplyr::select(cluster.nr, Sample, read.idx) -> read.idx.df
    
    left_join(fp.df, read.idx.df, by = c("Sample", "cluster.nr")) %>%
      unnest(read.idx) %>%
      unnest(start) %>%
      mutate(read.idx = as.factor(read.idx)) -> klee.df
    
    klee.df %>%
      mutate(biological.state = as.character(biological.state)) %>%
      mutate(biological.state = ifelse(biological.state %in% c("unrecognized"), "nucleosome", biological.state)) %>%
      mutate(biological.state = ifelse(biological.state %in% c("noise"), "accessible", biological.state)) %>%
      ggplot(aes(x = start, y = read.idx, fill=biological.state)) +
      geom_tile() +
      {if(is.null(partition.collapsing.dict)){facet_grid(cluster.nr~Sample, scales = "free_y", space = "free_y")}} +
      scale_x_continuous(limits = x.axis.breaks, breaks = x.axis.breaks, labels = format(x.axis.breaks, nsmall=1, big.mark=",")) +
      scale_fill_manual(breaks = c("TF", "accessible", "nucleosome"), values = c("dodgerblue4", "darkseagreen", "deepskyblue3")) +
      xlab(as.character(seqnames(footprint.charter.results$RegionOfInterest))) + 
      ylab(paste0(sum(distinct(fp.df, cluster.nr, read.count)$read.count), " molecules")) +
      theme_classic() +
      theme(text = element_text(size = 18), axis.title.y = element_text(vjust = 0), axis.title.x = element_text(vjust = 5),
            axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.border = element_rect(linewidth = 1, fill = "transparent"))
    
  }) -> klee.plots
  names(klee.plots) = current.sample
  
  return(klee.plots)  
  
}

#' @param data.type F1_bait_capture, F1_amplicon, WT_bait.capture, WT_amplicon, Rest_ko_amplicon, PRA_amplicon
#' @param deduplicate use only with PRA or amplicon data
#' @param plotting.TFBSs if NULL, defaults to isBound TFBSs in the 500bp window 
patch.single.site.plots = function(
    interpretable.master.table, rank, reference.genome = BSgenome.Mmusculus.UCSC.mm10, pool.replicates = FALSE, resize.size = 500, 
    k = NULL, partition.collapsing.dict = NULL, data.type = NULL, deduplicate = FALSE, plotting.TFBSs = NULL,
    remove.TFBS.labels = FALSE, reutrn.chromatin.influence.df = TRUE){
  
  Load.data(data.type = data.type, pool.replicates = pool.replicates)
  
  if(data.type %in% c("F1_bait.capture", "WT_bait.capture", "PRA_amplicon", "PRA_mutant_amplicon", "Sox2_kd_bait.capture", "Oct4_kd_bait.capture", "remodeller_bait.capture")){
    bait.capture.parameters.list$compute.pam.clustering["max.nr.reads"] = list(NULL)
    bait.capture.parameters.list$pam[["k"]] = ifelse(is.null(k), 8, k)
    bait.capture.parameters.list$output.params$plot.SW.heatmap = FALSE
    bait.capture.parameters.list$output.params$compute.silhouette = TRUE
    bait.capture.parameters.list$output.params$plot.binary.heatmap = TRUE
    params = bait.capture.parameters.list
  } else if (data.type %in% c("F1_amplicon", "WT_amplicon", "Rest_ko_amplicon", "Nrf1_kd_amplicon")){
    amplicon.parameters.list$compute.pam.clustering["max.nr.reads"] = list(NULL)
    amplicon.parameters.list$pam[["k"]] = ifelse(is.null(k), 10, k)
    amplicon.parameters.list$output.params$plot.SW.heatmap = FALSE
    amplicon.parameters.list$output.params$compute.silhouette = TRUE
    amplicon.parameters.list$output.params$plot.binary.heatmap = TRUE
    params = amplicon.parameters.list
  }
  
  current.sample = interpretable.master.table[rank,]$Sample
  current.tile = interpretable.master.table[rank,]$TFBS.cluster
  seqlevels(GenomicTiles) = seqlevels(TFBSs)
  RegionOfInterest = GenomicTiles[current.tile]
  RegionOfInterest_ext = IRanges::resize(RegionOfInterest, width = resize.size, fix = "center")
  
  if(data.type %in% c("F1_bait.capture", "F1_amplicon")){
    params$PlotSingleSiteSMF$PlottingSNPs = plyranges::filter_by_overlaps(SNPs[[ifelse(str_detect(current.sample, "^C"), "Cast", "Spret")]], RegionOfInterest_ext)
    if (length(params$PlotSingleSiteSMF$PlottingSNPs) == 0){
      params$PlotSingleSiteSMF$PlottingSNPs = NULL
    }
  } else if (data.type %in% c("WT_bait.capture", "WT_amplicon", "Rest_ko_amplicon", "Nrf1_kd_amplicon", "Sox2_kd_bait.capture", "Oct4_kd_bait.capture", "remodeller_bait.capture", "PRA_amplicon", "PRA_mutant_amplicon")) {
    params$PlotSingleSiteSMF$PlottingSNPs = NULL
  }
  
  if(data.type %in% c("F1_bait.capture", "F1_amplicon", "WT_bait.capture", "remodeller_bait.capture", "WT_amplicon", "Nrf1_kd_amplicon", "PRA_amplicon", "PRA_mutant_amplicon")){
    footprinting.type = "DE"
  } else if (data.type %in% c("Rest_ko_amplicon", "Sox2_kd_bait.capture", "Oct4_kd_bait.capture")){
    footprinting.type = "NO"
  }
  
  compute.pam.clustering.Wrapper(
    sampleSheet = sampleSheet, 
    Samples = MySamples[[footprinting.type]],
    RegionOfInterest = RegionOfInterest, 
    genome = reference.genome,
    TFBSs = TFBSs,
    params = params, 
    deduplicate = deduplicate,
    CytosinesToMask = CytosinesToMask,
    pool.samples = TRUE, # N.b: this affects which Cs get masked. FALSE give a species-specific mask, TRUE masks all Cs for all species. For GW I use TRUE, for single sites FALSE 
    SMF.plot = TRUE
  ) -> res
  
  if(data.type != "PRA_mutant_amplicon" & isTRUE(reutrn.chromatin.influence.df)){
    res %>% 
      create.master.table(
        list.of.files = ., single.site.coverage.thr = 20, cores = 1, 
        data.type = ifelse(data.type %in% c("WT_bait.capture", "WT_amplicon", "Rest_ko_amplicon", "Nrf1_kd_amplicon", "Sox2_kd_bait.capture", "Oct4_kd_bait.capture", "remodeller_bait.capture", "PRA_amplicon"), "PRA", "F1")
      ) %>%
      {if(data.type %in% c("F1")) filter(., Sample == current.sample) else .} %>%
      make.chromatin.influence.df(
        master.table = ., GenomicTiles = GenomicTiles, skip.chip = TRUE, 
        data.type = ifelse(data.type %in% c("WT_bait.capture", "WT_amplicon", "Rest_ko_amplicon", "Nrf1_kd_amplicon", "Sox2_kd_bait.capture", "Oct4_kd_bait.capture", "remodeller_bait.capture"), "WT", ifelse(data.type %in% c("PRA_amplicon"), "PRA", "F1"))
      ) -> chromatin.influence.df
  } else {
    chromatin.influence.df = NULL
  }
  
  if(data.type %in% c("F1_bait.capture", "F1_amplicon")){
    plotting.samples = sort(names(res$read.origin)[grepl(paste0("_", substr(current.sample,1,1)), gsub("_.*Barzaghi|_MethRate", "", names(res$read.origin)))], decreasing = TRUE)
  } else if (data.type %in% c("WT_bait.capture", "Rest_ko_amplicon", "Nrf1_kd_amplicon")){
    plotting.samples = names(res$read.origin)
  } else if (data.type %in% c("Sox2_kd_bait.capture", "Oct4_kd_bait.capture", "remodeller_bait.capture")){
    plotting.samples = rev(names(res$read.origin))
  } else if (data.type %in% c("WT_amplicon", "PRA_amplicon", "PRA_mutant_amplicon")){
    plotting.samples = grep(current.sample, names(res$read.origin), value = TRUE)
  }
  res$SingleSite.plot$data %<>%
    filter(sample %in% paste0(plotting.samples, "_MethRate")) %>% 
    mutate(sample = gsub("_.*Barzaghi|_MethRate", "", sample)) %>%
    mutate(sample = factor(sample, levels = plotting.samples))
  
  # bulk smf plot
  if(is.null(plotting.TFBSs)){plotting.TFBSs = plyranges::filter(plyranges::filter_by_overlaps(TFBSs, RegionOfInterest_ext), isBound)}
  polish.bulk.smf.plot(footprint.charter.results = res, data.type = data.type, current.sample = current.sample, RegionOfInterest = RegionOfInterest, plotting.TFBSs = plotting.TFBSs) -> res$SingleSite.plot
  if(isTRUE(remove.TFBS.labels)){res$SingleSite.plot$layers[[4]] = NULL}
  
  # SM stacks
  polish.sm.stacks(
    footprint.charter.results = res, current.sample = plotting.samples, 
    partition.collapsing.dict = partition.collapsing.dict
    ) -> SM.plots
  
  klee.pl = klee.plot(footprint.charter.results = res, current.sample = plotting.samples, partition.collapsing.dict = partition.collapsing.dict)
  
  if(data.type %in% c("F1_bait.capture", "F1_amplicon", "Nrf1_kd_amplicon", "Sox2_kd_bait.capture", "Oct4_kd_bait.capture", "remodeller_bait.capture")){
    res$SingleSite.plot + SM.plots + patchwork::plot_spacer() + klee.pl + 
      patchwork::plot_layout(
        ncol = 2, 
        guides = "collect",
        design = "
          A##
          BDC
          EDF
          ",
        heights = c(.5,.5),
        widths = c(.475,.05,.475)
      ) & theme(legend.box.just = "top") -> final.plot
    final.plot[[2]] = final.plot[[2]] + theme(axis.title.y = element_text(vjust = -8))
    final.plot[[5]] = final.plot[[5]] + theme(axis.title.y = element_text(vjust = -8))
    
  } else if (data.type %in% c("WT_bait.capture", "WT_amplicon", "Rest_ko_amplicon", "PRA_amplicon", "PRA_mutant_amplicon")){
    res$SingleSite.plot + SM.plots + patchwork::plot_spacer() + klee.pl + 
      patchwork::plot_layout(
        ncol = 2, 
        guides = "collect",
        design = "
          A##
          BCD
          ",
        heights = c(.5,.5),
        widths = c(.475,.05,.475)
      ) -> final.plot
    final.plot[[2]] = final.plot[[2]] + theme(axis.title.y = element_text(vjust = -8))
  }
  
  return(list(pl = final.plot, chromatin.influence.df = chromatin.influence.df))
  
}
