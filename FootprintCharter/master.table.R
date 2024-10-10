declare_dictionaries = function(){
  
  sequencing.batch.dict <<- list(
    batch.0 = c("SMF_MM_TKO"),
    batch.1 = c("CastIII", "CastIV", "SpretusIII", "SpretusIV", "Cast", "Spretus"),
    batch.2 = c("CITKO", "CIITKO", "SITKO", "SIITKO"),
    batch.3 = c("CIIITKO", "CIVTKO", "SIIITKO", "SIVTKO"),
    batch.4 = c("SVITKO", "SVIITKO"),
    batch.5 = c("SMF5", "SMF7", "SMF8", "SMF9", "SMFA", "SMFB", "SMFC", "SMFD"), # all AMP GB1
    batch.6 = c("S_5GB2", "S_7GB2", "4WTGB2", "5WTGB2", "S_1CS3", "S_2CS3", "1WTCS3", "2WTCS3", "3WTCS3"), # all AMP GB2,
    mixed = c("CTKO", "STKO") # when pooling samples
  )
  
  cell.line.dict <<- list(
    mESC.Spretus.F1 = c("SpretusIII", "SpretusIV", "Spretus"),
    mESC.Castaneus.F1 = c("CastIII", "CastIV", "Cast"),
    mESC.Spretus.F1.TKO = c("SITKO", "SIITKO", "SIIITKO", "SIVTKO", "SVITKO", "SVIITKO", "SMF5", "SMF7", "SMF8", "SMF9", "SMFA", "SMFB", "SMFC", "SMFD", "S_5GB2", "S_7GB2", "S_1CS3", "S_2CS3", "STKO"),
    mESC.Castaneus.F1.TKO = c("CITKO", "CIITKO", "CIIITKO", "CIVTKO", "CTKO"),
    mESC.WT.TKO = c("4WTGB2", "5WTGB2", "1WTCS3", "2WTCS3", "3WTCS3", "SMF_MM_TKO")
  )
  
  Footprinting.type.dict <<- list(
    DE = c("SITKO", "SIITKO", "SIIITKO", "SIVTKO", "SVITKO", "SVIITKO", "SMF5", "SMF7", "SMF8", "SMF9", "SMFA", "SMFB", "SMFC", "SMFD", "S_5GB2", "S_7GB2", "S_1CS3", "S_2CS3", 
           "CITKO", "CIITKO", "CIIITKO", "CIVTKO", 
           "4WTGB2", "5WTGB2", "1WTCS3", "2WTCS3", "3WTCS3",
           "CTKO", "STKO",
           "SMF_MM_TKO"),
    NO = c("CastIII", "CastIV", "SpretusIII", "SpretusIV", "Cast", "Spretus")
  )
  
  Sequencing.type.dict <<- list(
    bait.capture = c("CastIII", "CastIV", "SpretusIII", "SpretusIV", "SITKO", "SIITKO", "SIIITKO", "SIVTKO", "SVITKO", "SVIITKO", "CITKO", "CIITKO", "CIIITKO", "CIVTKO", "Cast", "Spretus", "CTKO", "STKO", "SMF_MM_TKO"),
    amplicon = c("SMF5", "SMF7", "SMF8", "SMF9", "SMFA", "SMFB", "SMFC", "SMFD", "S_5GB2", "S_7GB2", "4WTGB2", "5WTGB2", "S_1CS3", "S_2CS3", "1WTCS3", "2WTCS3", "3WTCS3")
  )
  
  clustering.batch.dict <<- list(
    batch.0 = c("SMF_MM_TKO"),
    batch.1 = c("SITKO", "SIITKO", "SIIITKO", "SIVTKO", "SVITKO", "SVIITKO", "CITKO", "CIITKO", "CIIITKO", "CIVTKO", "CTKO", "STKO"),
    batch.2 = c("CastIII", "CastIV", "SpretusIII", "SpretusIV", "Cast", "Spretus"),
    batch.3 = c("SMF5", "SMF7", "SMF8", "SMF9", "SMFA", "SMFB", "SMFC", "SMFD"),
    batch.4 = c("S_5GB2", "S_7GB2", "4WTGB2", "5WTGB2", "S_1CS3", "S_2CS3", "1WTCS3", "2WTCS3", "3WTCS3")
  )
  
}

#' @param list.of.unsupervised.clustering.results should be in the following form: list(results) | where results is coming from an Unsupervised.clustering.MultiSite.wrapper() call
#' @param partition.total.coverage DEPRECATED. Coverage threshold to pass (>=) to bring a partition down to statistical testing. Here by "total" we mean the summation of the Reference + Alternative alleles coverages.
#' @param single.site.coverage.thr defaults to 20. Minimum number of counts at EACH allele.
#' @param data.type F1 or PRA
create.master.table = function(list.of.files, 
                               partition.total.coverage = 0, 
                               single.site.coverage.thr = 20,
                               data.type = "F1",
                               cores = 1){
  
  # PREAMBLE
  if(is.list(list.of.files)){
    list.of.files$Footprinting.type = "DE"
    list.of.files$Sequencing.type = "bait.capture"
    list.of.files$clustering.batch = "batch.1"
    list.of.files$cluster.coverage.results.df = list.of.files$cluster.coverage.df
    list.of.files$unsupervised.footprints.detection.results.df = list.of.files$footprints.df
    list.of.files$unsupervised.footprints.detection.results.df %<>% mutate(TFBS.cluster = unique(list.of.files$cluster.coverage.df$TFBS.cluster))
    list.of.unsupervised.clustering.results = list(list.of.files)
  } else {
    source("/g/krebs/barzaghi/analyses/single.molecule.classification/utils/Source.all.unsupervised.functions.R")
    declare_dictionaries()
    list.of.files = unlist(list.of.files)
    list.of.unsupervised.clustering.results = lapply(seq_along(list.of.files), function(i){
      
      rslurm.res = readRDS(list.of.files[i])
      if(nrow(rslurm.res[[1]]$unsupervised.footprints.detection.results.df) == 0){return(NULL)}
      rslurm.res[[1]]$Footprinting.type = "DE"
      rslurm.res[[1]]$Sequencing.type = "bait.capture"
      rslurm.res[[1]]$clustering.batch = "batch.1"
      return(rslurm.res[[1]])
      
    })
    list.of.unsupervised.clustering.results = list.of.unsupervised.clustering.results[!unlist(lapply(list.of.unsupervised.clustering.results, is.null))]
  }
    
  if(partition.total.coverage>0){
    warning("parameter partition.total.coverage deprecated in favour of single.site.coverage.thr...")
    warning("Setting partition.total.coverage to 0")
    partition.total.coverage = 0
  }
  
  if(data.type == "F1"){
    
    Reduce(function(...) rbind(..., fill = TRUE), parallel::mclapply(seq_along(list.of.unsupervised.clustering.results), function(i){
      list.of.unsupervised.clustering.results[[i]]$cluster.coverage.results.df %>% 
        mutate(Footprinting.type = list.of.unsupervised.clustering.results[[i]]$Footprinting.type,
               Sequencing.type = list.of.unsupervised.clustering.results[[i]]$Sequencing.type,
               clustering.batch = list.of.unsupervised.clustering.results[[i]]$clustering.batch)
    }, mc.preschedule = TRUE, mc.cores = cores)) %>%
      filter(!is.na(cluster.nr) & read.count > 0) %>%
      mutate(cluster.nr = factor(cluster.nr, levels = seq(max(as.integer(cluster.nr), na.rm = TRUE)))) %>%
      separate(Sample, into = c("Allele", "Sample"), sep = "_", extra = "merge") %>%
      mutate(Allele = factor(Allele, levels = c("A", "R"))) %>%
      mutate(Sample = gsub("^.*_Barzaghi_|_dp_rm|Sample1|Sample2|^.*Kleinendorst_lane1|lane1|_$", "", Sample)) %>%
      mutate(Sample = gsub("^_|_$", "", gsub("DE|NO", "", Sample))) %>%
      distinct(cluster.nr, Allele, Sample, TFBS.cluster, .keep_all = TRUE) %>% #                                                                !!!!!!!!!!!! HERE YOU GOTTA FIX
      
      ##### Prepare data for statistical testing #####
    spread(Allele, read.count, drop = FALSE) %>% #                                                                                            <      N.b.: Here I already have 0 for the partitions that are not covered        >
      arrange(TFBS.cluster, Sample) %>% #                                                                                       <            NAs are only for the samples that are COMPLETELY uncovered           >
      filter(!is.na(A) & !is.na(R)) %>% #                                                                                       <     N.b.: at this line I'm collaterally discarding the mESC.WT.TKO Samples      >
      mutate(A = A+1, R = R+1) %>%
      group_by(Sample, TFBS.cluster) %>%
      mutate(tot.A = sum(A), tot.R = sum(R)) %>%
      filter(A+R >= partition.total.coverage) %>%
      filter(tot.A >= single.site.coverage.thr & tot.R >= single.site.coverage.thr) %>%
      ungroup() -> master.table
    
  } else if (data.type == "PRA"){
    
    Reduce(function(...) rbind(..., fill = TRUE), parallel::mclapply(seq_along(list.of.unsupervised.clustering.results), function(i){
      list.of.unsupervised.clustering.results[[i]]$cluster.coverage.results.df %>% 
        mutate(Footprinting.type = list.of.unsupervised.clustering.results[[i]]$Footprinting.type,
               Sequencing.type = list.of.unsupervised.clustering.results[[i]]$Sequencing.type,
               clustering.batch = list.of.unsupervised.clustering.results[[i]]$clustering.batch)
      }, mc.preschedule = TRUE, mc.cores = cores)) %>%
      filter(!is.na(cluster.nr) & read.count > 0) %>%
      mutate(cluster.nr = factor(cluster.nr, levels = seq(max(as.integer(cluster.nr), na.rm = TRUE)))) %>%
      mutate(Sample = gsub("^.*_Barzaghi_|_dp_rm|Sample1|Sample2|^.*Kleinendorst_lane1|lane1|_$", "", Sample)) %>%
      mutate(Sample = gsub("^_|_$", "", gsub("DE|NO", "", Sample))) %>%
      distinct(cluster.nr, Sample, TFBS.cluster, .keep_all = TRUE) %>% #                                                                !!!!!!!!!!!! HERE YOU GOTTA FIX

      ##### Prepare data for statistical testing #####
    arrange(TFBS.cluster, Sample) %>% #                                                                                       <            NAs are only for the samples that are COMPLETELY uncovered           >
      filter(!is.na(read.count)) %>% #                                                                                       <     N.b.: at this line I'm collaterally discarding the mESC.WT.TKO Samples      >
      mutate(read.count = read.count + 1) %>%
      group_by(Sample, TFBS.cluster) %>%
      mutate(tot.read.count = sum(read.count)) %>%
      filter(read.count >= partition.total.coverage) %>%
      filter(tot.read.count >= single.site.coverage.thr) %>%
      ungroup() -> master.table

  } else {
    stop("data.type not recognised")
  }
  
  if(nrow(master.table) == 0){
    
    return(master.table)
  
  } else {
    
    # Add biological.state interpretation to partitions
    Reduce(function(...) rbind(..., fill = TRUE), parallel::mclapply(seq_along(list.of.unsupervised.clustering.results), function(i){
      list.of.unsupervised.clustering.results[[i]]$unsupervised.footprints.detection.results.df %>% 
        mutate(Footprinting.type = list.of.unsupervised.clustering.results[[i]]$Footprinting.type, 
               Sequencing.type = list.of.unsupervised.clustering.results[[i]]$Sequencing.type, 
               clustering.batch = list.of.unsupervised.clustering.results[[i]]$clustering.batch)
    }, mc.preschedule = TRUE, mc.cores = cores)) %>%
      filter(!is.na(cluster.nr)) %>%
      mutate(cluster.nr = factor(cluster.nr, levels = seq(max(as.integer(cluster.nr), na.rm = TRUE)))) %>%
      dplyr::select(biological.state, cluster.nr, TFBS.cluster, Footprinting.type, Sequencing.type) %>%
      group_by(TFBS.cluster, cluster.nr, Footprinting.type, Sequencing.type) %>%
      summarise(biological.state = paste0(sort(unique(biological.state)), collapse = "+"),
                .groups = "drop") -> biological.state.table
    
    Reduce(function(...) rbind(..., fill = TRUE), parallel::mclapply(seq_along(list.of.unsupervised.clustering.results), function(i){
      list.of.unsupervised.clustering.results[[i]]$unsupervised.footprints.detection.results.df %>% 
        mutate(Footprinting.type = list.of.unsupervised.clustering.results[[i]]$Footprinting.type, 
               Sequencing.type = list.of.unsupervised.clustering.results[[i]]$Sequencing.type, 
               clustering.batch = list.of.unsupervised.clustering.results[[i]]$clustering.batch)
    }, mc.preschedule = TRUE, mc.cores = cores)) %>%
      filter(!is.na(cluster.nr)) %>%
      mutate(cluster.nr = factor(cluster.nr, levels = seq(max(as.integer(cluster.nr), na.rm = TRUE)))) %>%
      filter(biological.state == "TF") %>%
      mutate(footprint.center = as.integer(format(round(rowMeans(.[,c("start","end")])), scientific=FALSE)),
             footprint.width = end - start + 1) %>%
      dplyr::select(TF, TF.name, 
             footprint.center, footprint.width,
             footprint.idx, cluster.nr, TFBS.cluster, Footprinting.type, Sequencing.type) %>%
      rowwise() %>%
      mutate(TF = paste0(TF,".",footprint.idx, collapse = "+")) %>%
      mutate(TF.name = paste0(TF.name, collapse = "+")) %>%
      mutate(footprint.center = paste0(rep(footprint.center, str_count(TF, "\\+")+1), collapse = "+")) %>%
      mutate(footprint.width = paste0(rep(footprint.width, str_count(TF, "\\+")+1), collapse = "+")) %>%
      ungroup() %>%
      distinct() %>%
      group_by(TFBS.cluster, cluster.nr, Footprinting.type, Sequencing.type) %>%
      summarise(
        TF = paste0(TF, collapse = "+"), 
        TF.name = paste0(TF.name, collapse = "+"),
        footprint.center = paste0(footprint.center, collapse = "+"),
        footprint.width = paste0(footprint.width, collapse = "+"),
        .groups = "drop") -> TF.table
    
    Reduce(function(...) rbind(..., fill = TRUE), parallel::mclapply(seq_along(list.of.unsupervised.clustering.results), function(i){
      list.of.unsupervised.clustering.results[[i]]$unsupervised.footprints.detection.results.df %>% 
        mutate(Footprinting.type = list.of.unsupervised.clustering.results[[i]]$Footprinting.type, 
               Sequencing.type = list.of.unsupervised.clustering.results[[i]]$Sequencing.type, 
               clustering.batch = list.of.unsupervised.clustering.results[[i]]$clustering.batch)
    }, mc.preschedule = TRUE, mc.cores = cores)) %>%
      filter(!is.na(cluster.nr)) %>%
      mutate(cluster.nr = factor(cluster.nr, levels = seq(max(as.integer(cluster.nr), na.rm = TRUE)))) %>%
      filter(biological.state == "nucleosome") %>%
      mutate(nucleosome_footprint.center = as.integer(format(round(rowMeans(.[,c("start","end")])), scientific=FALSE)),
             nucleosome_footprint.width = end - start + 1) %>%
      dplyr::select(nucleosome_footprint.center, nucleosome_footprint.width,
             footprint.idx, cluster.nr, TFBS.cluster, Footprinting.type, Sequencing.type) %>%
      rowwise() %>%
      mutate(nucleosome_footprint.center = paste0(nucleosome_footprint.center, collapse = "+")) %>%
      mutate(nucleosome_footprint.width = paste0(nucleosome_footprint.width, collapse = "+")) %>%
      ungroup() %>%
      distinct() %>%
      group_by(TFBS.cluster, cluster.nr, Footprinting.type, Sequencing.type) %>%
      summarise(
        nucleosome_footprint.center = paste0(nucleosome_footprint.center, collapse = "+"),
        nucleosome_footprint.width = paste0(nucleosome_footprint.width, collapse = "+"),
        .groups = "drop") -> nucleosome.table
    
    Reduce(function(...) rbind(..., fill = TRUE), parallel::mclapply(seq_along(list.of.unsupervised.clustering.results), function(i){
      list.of.unsupervised.clustering.results[[i]]$unsupervised.footprints.detection.results.df %>% 
        mutate(Footprinting.type = list.of.unsupervised.clustering.results[[i]]$Footprinting.type, 
               Sequencing.type = list.of.unsupervised.clustering.results[[i]]$Sequencing.type, 
               clustering.batch = list.of.unsupervised.clustering.results[[i]]$clustering.batch)
    }, mc.preschedule = TRUE, mc.cores = cores)) %>%
      filter(!is.na(cluster.nr)) %>%
      mutate(cluster.nr = factor(cluster.nr, levels = seq(max(as.integer(cluster.nr), na.rm = TRUE)))) %>%
      filter(biological.state == "unrecognized") %>%
      mutate(unrecognized_footprint.center = as.integer(format(round(rowMeans(.[,c("start","end")])), scientific=FALSE)),
             unrecognized_footprint.width = end - start + 1) %>%
      dplyr::select(unrecognized_footprint.center, unrecognized_footprint.width,
             footprint.idx, cluster.nr, TFBS.cluster, Footprinting.type, Sequencing.type) %>%
      rowwise() %>%
      mutate(unrecognized_footprint.center = paste0(unrecognized_footprint.center, collapse = "+")) %>%
      mutate(unrecognized_footprint.width = paste0(unrecognized_footprint.width, collapse = "+")) %>%
      ungroup() %>%
      distinct() %>%
      group_by(TFBS.cluster, cluster.nr, Footprinting.type, Sequencing.type) %>%
      summarise(
        unrecognized_footprint.center = paste0(unrecognized_footprint.center, collapse = "+"),
        unrecognized_footprint.width = paste0(unrecognized_footprint.width, collapse = "+"),
        .groups = "drop") -> unrecognized.table
    
    Reduce(function(...) rbind(..., fill = TRUE), parallel::mclapply(seq_along(list.of.unsupervised.clustering.results), function(i){
      list.of.unsupervised.clustering.results[[i]]$unsupervised.footprints.detection.results.df %>% 
        mutate(Footprinting.type = list.of.unsupervised.clustering.results[[i]]$Footprinting.type, 
               Sequencing.type = list.of.unsupervised.clustering.results[[i]]$Sequencing.type, 
               clustering.batch = list.of.unsupervised.clustering.results[[i]]$clustering.batch)
    }, mc.preschedule = TRUE, mc.cores = cores)) %>%
      filter(!is.na(cluster.nr)) %>%
      mutate(cluster.nr = factor(cluster.nr, levels = seq(max(as.integer(cluster.nr), na.rm = TRUE)))) %>%
      filter(biological.state == "noise") %>%
      mutate(noise_footprint.center = as.integer(format(round(rowMeans(.[,c("start","end")])), scientific=FALSE)),
             noise_footprint.width = end - start + 1) %>%
      dplyr::select(noise_footprint.center, noise_footprint.width,
             footprint.idx, cluster.nr, TFBS.cluster, Footprinting.type, Sequencing.type) %>%
      rowwise() %>%
      mutate(noise_footprint.center = paste0(noise_footprint.center, collapse = "+")) %>%
      mutate(noise_footprint.width = paste0(noise_footprint.width, collapse = "+")) %>%
      ungroup() %>%
      distinct() %>%
      group_by(TFBS.cluster, cluster.nr, Footprinting.type, Sequencing.type) %>%
      summarise(
        noise_footprint.center = paste0(noise_footprint.center, collapse = "+"),
        noise_footprint.width = paste0(noise_footprint.width, collapse = "+"),
        .groups = "drop") -> noise.table
    
    Reduce(function(...) rbind(..., fill = TRUE), parallel::mclapply(seq_along(list.of.unsupervised.clustering.results), function(i){
      list.of.unsupervised.clustering.results[[i]]$unsupervised.footprints.detection.results.df %>% 
        mutate(Footprinting.type = list.of.unsupervised.clustering.results[[i]]$Footprinting.type, 
               Sequencing.type = list.of.unsupervised.clustering.results[[i]]$Sequencing.type, 
               clustering.batch = list.of.unsupervised.clustering.results[[i]]$clustering.batch)
    }, mc.preschedule = TRUE, mc.cores = cores)) %>%
      filter(!is.na(cluster.nr)) %>%
      mutate(cluster.nr = factor(cluster.nr, levels = seq(max(as.integer(cluster.nr), na.rm = TRUE)))) %>%
      filter(biological.state == "accessible") %>%
      mutate(accessible_footprint.center = as.integer(format(round(rowMeans(.[,c("start","end")])), scientific=FALSE)),
             accessible_footprint.width = end - start + 1) %>%
      dplyr::select(accessible_footprint.center, accessible_footprint.width,
             footprint.idx, cluster.nr, TFBS.cluster, Footprinting.type, Sequencing.type) %>%
      rowwise() %>%
      mutate(accessible_footprint.center = paste0(accessible_footprint.center, collapse = "+")) %>%
      mutate(accessible_footprint.width = paste0(accessible_footprint.width, collapse = "+")) %>%
      ungroup() %>%
      distinct() %>%
      group_by(TFBS.cluster, cluster.nr, Footprinting.type, Sequencing.type) %>%
      summarise(
        accessible_footprint.center = paste0(accessible_footprint.center, collapse = "+"),
        accessible_footprint.width = paste0(accessible_footprint.width, collapse = "+"),
        .groups = "drop") -> accessibility.table
    
    # Join everything
    full_join(TF.table, nucleosome.table, by = c("TFBS.cluster", "cluster.nr", "Footprinting.type", "Sequencing.type")) -> footprints.table_tf.nuc
    full_join(footprints.table_tf.nuc, unrecognized.table, by = c("TFBS.cluster", "cluster.nr", "Footprinting.type", "Sequencing.type")) -> footprints.table_tf.nuc.unr
    full_join(footprints.table_tf.nuc.unr, noise.table, by = c("TFBS.cluster", "cluster.nr", "Footprinting.type", "Sequencing.type")) -> footprints.table_tf.nuc.unr.noi.
    full_join(footprints.table_tf.nuc.unr.noi., accessibility.table, by = c("TFBS.cluster", "cluster.nr", "Footprinting.type", "Sequencing.type")) -> footprints.table_tf.nuc.unr.noi.acc
    full_join(biological.state.table, footprints.table_tf.nuc.unr.noi.acc, by = c("TFBS.cluster", "cluster.nr", "Footprinting.type", "Sequencing.type")) -> biological.interpretation.table
    full_join(master.table, biological.interpretation.table, by = c("TFBS.cluster", "cluster.nr", "Footprinting.type", "Sequencing.type")) -> interpretable.master.table
    
    interpretable.master.table = ungroup(interpretable.master.table)
    
    return(interpretable.master.table)
    
  }
    
}
