annotated_tiles_without_footprints = function(master.table, GenomicTiles, TFBSs, data.type = "F1"){

  master.table %>%
    mutate(TF.name = ifelse(nchar(gsub("NA|\\+", "", TF.name)) == 0, NA, TF.name)) %>%
    group_by(Sample, TFBS.cluster) %>%
    filter(all(is.na(TF.name))) %>% # filter out the whole Tile if there is a detected footprint. Add  if wanting to retain unannotated footprints
    ungroup() -> tmp

  plyranges::find_overlaps(
    TFBSs, # %>% plyranges::filter(isBound),
    GenomicTiles[unique(tmp$TFBS.cluster)] %>% plyranges::mutate(TFBS.cluster = names(.))
    ) %>%
    plyranges::select(TF, absolute.idx, TFBS.cluster) %>%
    as_tibble() %>%
    dplyr::select(TF, absolute.idx, TFBS.cluster) %>%
    mutate(TF = paste0(TF, ".0")) %>% # placeholder to match master.table structure
    group_by(TFBS.cluster) %>%
    summarise(TF = paste(TF, collapse = "+"), TF.name = paste(absolute.idx, collapse = "+")) -> annotation

  left_join(
    tmp %>% select(-TF, -TF.name),
    annotation,
    by = "TFBS.cluster"
  ) %>%
    mutate(R = NA, A = NA, footprint.center = NA, footprint.width = NA) -> tmp.2

  if(data.type %in% c("PRA", "WT")){
    tmp.2 %<>%
      select(-A, -R) %>%
      mutate(read.count = NA)
  }

  master.table %>%
    group_by(Sample, TFBS.cluster) %>%
    filter(any(!is.na(TF.name)) & any(nchar(gsub("NA|\\+", "", TF.name)) > 0)) %>% # filter out the whole Tile if there is a detected footprint. Add  if wanting to retain unannotated footprints
    ungroup()  -> tmp.3

  rbind(tmp.2, tmp.3) -> return.df

  return(return.df)

}

#' @param TFBSs has to have an isBound boolean variable encoding ChIP{seq,nexus,exo} bound state
#' @param normalise.by.IC whether to normalise PWM scores and deltas by the information content of the relative PWM
add.delta.pwm = function(master.table, TFBSs, normalise.by.IC = FALSE, data.type = "F1"){

  warning("isBound coulumn used to filter for truly bound sites")
  warning("N.b.: the !is.na(delta.pwm) filter (given the plyranges::filter(isBound) filter above) imposes the we are looking only at ChIPseq bound sites")
  warning("[UPDATE 02.03.23] we're keeping also the non footprintables")

  TFBSs_covered = unique(unlist(str_split(unique(master.table$TF.name), "\\+")))

  delta.pwm.df = TFBSs %>%
    plyranges::filter(absolute.idx %in% TFBSs_covered) %>% # isBound, 
    as_tibble() %>%
    dplyr::select(TF, absolute.idx, BL6.absScore, Cast.delta.pwm, Spret.delta.pwm) %>% # N.b.: these deltas are A-R (I checked)
    dplyr::rename(CTKO = "Cast.delta.pwm", STKO = "Spret.delta.pwm") %>%
    dplyr::rename(TF.name = "absolute.idx") %>%
    gather(Sample, delta.pwm, CTKO, STKO) %>%
    mutate(delta.pwm = as.numeric(delta.pwm))

  if (data.type == "F1") {
    
    # needed to ducktape normalise.by.IC
    delta.pwm.df %<>%
      select(-TF)

    left_join(
        master.table %>%
          mutate(original.sample = Sample, Sample = gsub("I|V", "", Sample)) %>% # deals with the different samples names when processing unpooled replicates
          mutate(TF = str_split(TF, "\\+"), TF.name = str_split(TF.name, "\\+"), footprint.center = str_split(footprint.center, "\\+"), footprint.width = str_split(footprint.width, "\\+")) %>%
          unnest(c(TF, TF.name, footprint.center, footprint.width)) %>%
          mutate(footprint.center = as.integer(footprint.center), footprint.width = as.integer(footprint.width)) %>%
          filter(TF.name != "NA") %>%
          mutate(footprint.start = ifelse(is.na(footprint.width) & is.na(footprint.center), start(TFBSs[TF.name]), round(footprint.center - footprint.width/2)),
                 footprint.end = ifelse(is.na(footprint.width) & is.na(footprint.center), end(TFBSs[TF.name]), round(footprint.center + footprint.width/2))) %>%
          group_by(original.sample, Sample, TFBS.cluster, CGI, cluster.nr, tot.A, tot.R, TF, TF.name) %>%
          summarise(A = sum(A), R = sum(R), cluster.nr.list = paste0(cluster.nr, collapse = "+"),
                    footprint.width = max(footprint.end, na.rm = TRUE) - min(footprint.start, na.rm = TRUE), # ifelse(all(is.na(footprint.width)), NA,  max(footprint.end, na.rm = TRUE) - min(footprint.start, na.rm = TRUE)),
                    footprint.center = min(footprint.start, na.rm = TRUE) + footprint.width/2, #ifelse(all(is.na(footprint.center)), NA, min(footprint.start, na.rm = TRUE) + footprint.width/2), # N.b.: here we melt footprints into the largest possible one
                    .groups = "drop") %>%
          mutate(cluster.nr.list = ifelse(is.na(A) & is.na(R), NA, cluster.nr.list)),
        delta.pwm.df,
        by = c("Sample", "TF.name")
      ) %>%
        mutate(Sample = original.sample) %>% # deals with the different samples names when processing unpooled replicates
        select(-original.sample) %>% # deals with the different samples names when processing unpooled replicates
        filter(!is.na(delta.pwm)) -> master.table.delta.pwm

  } else if (data.type %in% c("PRA", "WT")){

    left_join(
      master.table %>%
        mutate(CGI = NA) %>%
        mutate(original.sample = Sample, Sample = gsub("I|V", "", Sample)) %>% # deals with the different samples names when processing unpooled replicates
        mutate(TF = str_split(TF, "\\+"), TF.name = str_split(TF.name, "\\+"), footprint.center = str_split(footprint.center, "\\+"), footprint.width = str_split(footprint.width, "\\+")) %>%
        unnest(c(TF, TF.name, footprint.center, footprint.width)) %>%
        mutate(footprint.center = as.integer(footprint.center), footprint.width = as.integer(footprint.width)) %>%
        filter(TF.name != "NA") %>%
        mutate(footprint.start = ifelse(is.na(footprint.width) & is.na(footprint.center), start(TFBSs[TF.name]), round(footprint.center - footprint.width/2)),
               footprint.end = ifelse(is.na(footprint.width) & is.na(footprint.center), end(TFBSs[TF.name]), round(footprint.center + footprint.width/2))) %>%
        group_by(original.sample, Sample, TFBS.cluster, CGI, cluster.nr, tot.read.count, TF, TF.name) %>%
        summarise(read.count = sum(read.count), cluster.nr.list = paste0(cluster.nr, collapse = "+"),
                  footprint.width = max(footprint.end, na.rm = TRUE) - min(footprint.start, na.rm = TRUE), # ifelse(all(is.na(footprint.width)), NA,  max(footprint.end, na.rm = TRUE) - min(footprint.start, na.rm = TRUE)),
                  footprint.center = min(footprint.start, na.rm = TRUE) + footprint.width/2, #ifelse(all(is.na(footprint.center)), NA, min(footprint.start, na.rm = TRUE) + footprint.width/2), # N.b.: here we melt footprints into the largest possible one
                  .groups = "drop") %>%
        mutate(cluster.nr.list = ifelse(is.na(read.count), NA, cluster.nr.list)) %>%
        mutate(Sample = original.sample) %>% # deals with the different samples names when processing unpooled replicates
        select(-original.sample), # deals with the different samples names when processing unpooled replicates
      delta.pwm.df %>%
        select(TF.name, BL6.absScore) %>%
        distinct(),
      by = "TF.name"
    ) %>%
      mutate(delta.pwm = 0) -> master.table.delta.pwm

  }

  return(master.table.delta.pwm)

}

filter.bound.TFBSs = function(x, TFBSs){
  
  return(
    filter(x, TF.name %in% plyranges::filter(TFBSs, isBound)$absolute.idx)
  )
  
}

#' if two footprints are annotated with same TF & TF.name, aggregate their count: i.e. I don't believe these sub-footprints to be biologically relevant
aggregate.footprints.by.TF.annotation = function(master.table, TFBSs, data.type = "F1"){

  if(data.type == "F1"){

    master.table %>%
      separate(TF, into = c("TF", "fp.idx"), sep = "\\.") %>%
      mutate(footprint.start = ifelse(is.na(footprint.width) & is.na(footprint.center), start(TFBSs[TF.name]), round(footprint.center - footprint.width/2)),
             footprint.end = ifelse(is.na(footprint.width) & is.na(footprint.center), end(TFBSs[TF.name]), round(footprint.center + footprint.width/2))) %>%
      # Aggregate footprints by TF annotation
      # group_by(Sample, TFBS.cluster, CGI, cluster.nr, tot.A, tot.R, TF, TF.name, score, delta.pwm, flanking.SNP) %>%
      group_by(Sample, TFBS.cluster, CGI, cluster.nr, tot.A, tot.R, TF, TF.name, BL6.absScore, delta.pwm) %>%
      summarise(fp.idx = min(fp.idx),
                A = unique(A), R=unique(R), cluster.nr.list = paste0(cluster.nr, collapse = "+"),
                footprint.width = max(footprint.end, na.rm = TRUE) - min(footprint.start, na.rm = TRUE), # ifelse(all(is.na(footprint.width)), NA,  max(footprint.end, na.rm = TRUE) - min(footprint.start, na.rm = TRUE)),
                footprint.center = min(footprint.start, na.rm = TRUE) + footprint.width/2, # ifelse(all(is.na(footprint.center)), NA, min(footprint.start, na.rm = TRUE) + footprint.width/2), # N.b.: here we melt footprints into the largest possible one
                .groups = "drop") %>%
      # Aggregate counts across clusters
      mutate(footprint.start = ifelse(is.na(footprint.width) & is.na(footprint.center), start(TFBSs[TF.name]), round(footprint.center - footprint.width/2)),
             footprint.end = ifelse(is.na(footprint.width) & is.na(footprint.center), end(TFBSs[TF.name]), round(footprint.center + footprint.width/2))) %>%
      # group_by(Sample, TFBS.cluster, CGI, tot.A, tot.R, TF, TF.name, score, delta.pwm, flanking.SNP) %>%
      group_by(Sample, TFBS.cluster, CGI, tot.A, tot.R, TF, TF.name, BL6.absScore, delta.pwm) %>%
      summarise(fp.idx = min(fp.idx),
                A = sum(A), R=sum(R), cluster.nr.list = paste0(cluster.nr, collapse = "+"),
                footprint.width = max(footprint.end, na.rm = TRUE) - min(footprint.start, na.rm = TRUE), # ifelse(all(is.na(footprint.width)), NA,  max(footprint.end, na.rm = TRUE) - min(footprint.start, na.rm = TRUE)),
                footprint.center = min(footprint.start, na.rm = TRUE) + footprint.width/2, # ifelse(all(is.na(footprint.center)), NA, min(footprint.start, na.rm = TRUE) + footprint.width/2), # N.b.: here we melt footprints into the largest possible one
                .groups = "drop") %>%
      mutate(cluster.nr.list = ifelse(is.na(A) & is.na(R), NA, cluster.nr.list)) %>%
      unite("TF", c(TF, fp.idx), sep = ".")

  } else if (data.type %in% c("PRA", "WT")){

    master.table %>%
      separate(TF, into = c("TF", "fp.idx"), sep = "\\.") %>%
      mutate(footprint.start = ifelse(is.na(footprint.width) & is.na(footprint.center), start(TFBSs[TF.name]), round(footprint.center - footprint.width/2)),
             footprint.end = ifelse(is.na(footprint.width) & is.na(footprint.center), end(TFBSs[TF.name]), round(footprint.center + footprint.width/2))) %>%
      # Aggregate footprints by TF annotation
      # group_by(Sample, TFBS.cluster, CGI, cluster.nr, tot.A, tot.R, TF, TF.name, score, delta.pwm, flanking.SNP) %>%
      group_by(Sample, TFBS.cluster, CGI, cluster.nr, tot.read.count, TF, TF.name, BL6.absScore, delta.pwm) %>%
      summarise(fp.idx = min(fp.idx),
                read.count = unique(read.count), cluster.nr.list = paste0(cluster.nr, collapse = "+"),
                footprint.width = max(footprint.end, na.rm = TRUE) - min(footprint.start, na.rm = TRUE), # ifelse(all(is.na(footprint.width)), NA,  max(footprint.end, na.rm = TRUE) - min(footprint.start, na.rm = TRUE)),
                footprint.center = min(footprint.start, na.rm = TRUE) + footprint.width/2, # ifelse(all(is.na(footprint.center)), NA, min(footprint.start, na.rm = TRUE) + footprint.width/2), # N.b.: here we melt footprints into the largest possible one
                .groups = "drop") %>%
      # Aggregate counts across clusters
      mutate(footprint.start = ifelse(is.na(footprint.width) & is.na(footprint.center), start(TFBSs[TF.name]), round(footprint.center - footprint.width/2)),
             footprint.end = ifelse(is.na(footprint.width) & is.na(footprint.center), end(TFBSs[TF.name]), round(footprint.center + footprint.width/2))) %>%
      # group_by(Sample, TFBS.cluster, CGI, tot.A, tot.R, TF, TF.name, score, delta.pwm, flanking.SNP) %>%
      group_by(Sample, TFBS.cluster, CGI, tot.read.count, TF, TF.name, BL6.absScore, delta.pwm) %>%
      summarise(fp.idx = min(fp.idx),
                read.count=sum(read.count), cluster.nr.list = paste0(cluster.nr, collapse = "+"),
                footprint.width = max(footprint.end, na.rm = TRUE) - min(footprint.start, na.rm = TRUE), # ifelse(all(is.na(footprint.width)), NA,  max(footprint.end, na.rm = TRUE) - min(footprint.start, na.rm = TRUE)),
                footprint.center = min(footprint.start, na.rm = TRUE) + footprint.width/2, # ifelse(all(is.na(footprint.center)), NA, min(footprint.start, na.rm = TRUE) + footprint.width/2), # N.b.: here we melt footprints into the largest possible one
                .groups = "drop") %>%
      mutate(cluster.nr.list = ifelse(is.na(read.count), NA, cluster.nr.list)) %>%
      unite("TF", c(TF, fp.idx), sep = ".")

  }


}

filter.TFs.for.binding.frequency = function(master.table, THR, x.type = "singleTF", data.type = "F1"){

  if(x.type == "singleTF"){

    if(data.type == "F1"){

      master.table %>%
        mutate(
          A = ifelse(A/tot.A + R/tot.R <= THR, NA, A),
          R = ifelse(A/tot.A + R/tot.R <= THR, NA, R)#,
          # footprint.width = ifelse(A/tot.A <= THR & R/tot.R <= THR, NA, footprint.width),
          # footprint.center = ifelse(A/tot.A <= THR & R/tot.R <= THR, NA, footprint.center),
        ) %>%
        mutate(cluster.nr.list = ifelse(is.na(A) & is.na(R), NA, cluster.nr.list))

    } else if (data.type %in% c("PRA", "WT")){

      master.table %>%
        mutate(
          read.count = ifelse(read.count/tot.read.count <= THR, NA, read.count)
        ) %>%
        mutate(cluster.nr.list = ifelse(is.na(read.count), NA, cluster.nr.list))

    }


  } else if (x.type == "TFpairs"){

    if(data.type == "F1"){

      master.table %>%
        mutate(
          A = ifelse(!(R/tot.R + A/tot.A > THR | is.na(footprint.width)) & (cobinder.R/cobinder.tot.R + cobinder.A/cobinder.tot.A > THR), NA, A),
          R = ifelse(!(R/tot.R + A/tot.A > THR | is.na(footprint.width)) & (cobinder.R/cobinder.tot.R + cobinder.A/cobinder.tot.A > THR), NA, R),
          cobinder.R = ifelse(!(R/tot.R + A/tot.A > THR | is.na(footprint.width)) & (cobinder.R/cobinder.tot.R + cobinder.A/cobinder.tot.A > THR), NA, cobinder.R),
          cobinder.A = ifelse(!(R/tot.R + A/tot.A > THR | is.na(footprint.width)) & (cobinder.R/cobinder.tot.R + cobinder.A/cobinder.tot.A > THR), NA, cobinder.A)#,
          # footprint.width = ifelse(!(R/tot.R > THR | A/tot.A > THR | is.na(footprint.width)) & (cobinder.R/cobinder.tot.R > THR | cobinder.A/cobinder.tot.A > THR), NA, footprint.width),
          # footprint.center = ifelse(!(R/tot.R > THR | A/tot.A > THR | is.na(footprint.width)) & (cobinder.R/cobinder.tot.R > THR | cobinder.A/cobinder.tot.A > THR), NA, footprint.center),
          # cobinder.footprint.width = ifelse(!(R/tot.R > THR | A/tot.A > THR | is.na(footprint.width)) & (cobinder.R/cobinder.tot.R > THR | cobinder.A/cobinder.tot.A > THR), NA, cobinder.footprint.width),
          # cobinder.footprint.center = ifelse(!(R/tot.R > THR | A/tot.A > THR | is.na(footprint.width)) & (cobinder.R/cobinder.tot.R > THR | cobinder.A/cobinder.tot.A > THR), NA, cobinder.footprint.center)
        ) %>%
        mutate(cluster.nr.list = ifelse(is.na(A) & is.na(R), NA, cluster.nr.list))

    } else if (data.type %in% c("PRA", "WT")){

      stop("YOU NEED TO IMPLEMENT THIS")

    }



  }

}

#' @param x data.frame to add nucleosomal occupancy info to
#' @param master.table unprocessed and unfiltered master table
#' @param window.width if NULL (default) width will be the TF motif/footprint. Otherwise pass integer
add.nucleosome.occupancy = function(x, master.table, data.type = "F1", window.width = NULL){

  if(data.type == "F1"){

    master.table %>%
      filter(Sample %in% unique(x$Sample), TFBS.cluster %in% unique(x$TFBS.cluster), str_detect(biological.state, "nucleosome")) %>%
      dplyr::select(cluster.nr, Sample, TFBS.cluster, tot.A, tot.R, A, R, nucleosome_footprint.center, nucleosome_footprint.width) %>%
      mutate(center = str_split(nucleosome_footprint.center, "\\+"), width = str_split(nucleosome_footprint.width, "\\+")) %>%
      dplyr::select(-nucleosome_footprint.center, -nucleosome_footprint.width) %>%
      unnest(c(center, width)) %>%
      mutate(center = as.integer(center), width = as.integer(width)) %>%
      mutate(start = center - width/2, end = center + width/2) %>%
      dplyr::select(-center, -width) %>%
      unite("seqnames", c("Sample", "TFBS.cluster"), sep = ";") %>%
      plyranges::as_granges() -> nuc.granges
    
    master.table %>%
      filter(Sample %in% unique(x$Sample), TFBS.cluster %in% unique(x$TFBS.cluster), str_detect(biological.state, "unrecognized")) %>%
      dplyr::select(cluster.nr, Sample, TFBS.cluster, tot.A, tot.R, A, R, unrecognized_footprint.center, unrecognized_footprint.width) %>%
      mutate(center = str_split(unrecognized_footprint.center, "\\+"), width = str_split(unrecognized_footprint.width, "\\+")) %>%
      dplyr::select(-unrecognized_footprint.center, -unrecognized_footprint.width) %>%
      unnest(c(center, width)) %>%
      mutate(center = as.integer(center), width = as.integer(width)) %>%
      mutate(start = center - width/2, end = center + width/2) %>%
      dplyr::select(-center, -width) %>%
      unite("seqnames", c("Sample", "TFBS.cluster"), sep = ";") %>%
      plyranges::as_granges() -> unr.granges

    master.table %>%
      filter(Sample %in% unique(x$Sample), TFBS.cluster %in% unique(x$TFBS.cluster), str_detect(biological.state, "accessible|noise")) %>%
      unite("accessible_footprint.center", c("accessible_footprint.center", "footprint.center", "noise_footprint.center"), sep = "+") %>% # Warning: unsorted | this adds the TF footprints counts
      unite("accessible_footprint.width", c("accessible_footprint.width", "footprint.width", "noise_footprint.width"), sep = "+") %>% # Warning: unsorted | this adds the TF footprints counts
      # mutate(accessible_footprint.center = gsub("^NA\\+|\\+NA$", "", accessible_footprint.center), accessible_footprint.width = gsub("^NA\\+|\\+NA$", "", accessible_footprint.width)) %>%
      # mutate(accessible_footprint.center = gsub("\\+NA", "", accessible_footprint.center), accessible_footprint.width = gsub("\\+NA", "", accessible_footprint.width)) %>%
      dplyr::select(cluster.nr, Sample, TFBS.cluster, tot.A, tot.R, A, R, accessible_footprint.center, accessible_footprint.width) %>%
      mutate(center = str_split(accessible_footprint.center, "\\+"), width = str_split(accessible_footprint.width, "\\+")) %>%
      dplyr::select(-accessible_footprint.center, -accessible_footprint.width) %>%
      unnest(c(center, width)) %>%
      filter(center!="NA" & width!="NA") %>%
      mutate(center = as.integer(center), width = as.integer(width)) %>%
      mutate(start = center - width/2, end = center + width/2) %>%
      dplyr::select(-center, -width) %>%
      unite("seqnames", c("Sample", "TFBS.cluster"), sep = ";") %>%
      plyranges::as_granges() -> acc.granges

  } else if (data.type %in% c("PRA", "WT")){

    master.table %>%
      filter(Sample %in% unique(x$Sample), TFBS.cluster %in% unique(x$TFBS.cluster), str_detect(biological.state, "nucleosome")) %>%
      dplyr::select(cluster.nr, Sample, TFBS.cluster, tot.read.count, read.count, nucleosome_footprint.center, nucleosome_footprint.width) %>%
      mutate(center = str_split(nucleosome_footprint.center, "\\+"), width = str_split(nucleosome_footprint.width, "\\+")) %>%
      dplyr::select(-nucleosome_footprint.center, -nucleosome_footprint.width) %>%
      unnest(c(center, width)) %>%
      mutate(center = as.integer(center), width = as.integer(width)) %>%
      mutate(start = center - width/2, end = center + width/2) %>%
      dplyr::select(-center, -width) %>%
      unite("seqnames", c("Sample", "TFBS.cluster"), sep = ";") %>%
      plyranges::as_granges() -> nuc.granges
    
    master.table %>%
      filter(Sample %in% unique(x$Sample), TFBS.cluster %in% unique(x$TFBS.cluster), str_detect(biological.state, "unrecognized")) %>%
      dplyr::select(cluster.nr, Sample, TFBS.cluster, tot.read.count, read.count, unrecognized_footprint.center, unrecognized_footprint.width) %>%
      mutate(center = str_split(unrecognized_footprint.center, "\\+"), width = str_split(unrecognized_footprint.width, "\\+")) %>%
      dplyr::select(-unrecognized_footprint.center, -unrecognized_footprint.width) %>%
      unnest(c(center, width)) %>%
      mutate(center = as.integer(center), width = as.integer(width)) %>%
      mutate(start = center - width/2, end = center + width/2) %>%
      dplyr::select(-center, -width) %>%
      unite("seqnames", c("Sample", "TFBS.cluster"), sep = ";") %>%
      plyranges::as_granges() -> unr.granges
    
    master.table %>%
      filter(Sample %in% unique(x$Sample), TFBS.cluster %in% unique(x$TFBS.cluster), str_detect(biological.state, "accessible|noise")) %>%
      unite("accessible_footprint.center", c("accessible_footprint.center", "footprint.center", "noise_footprint.center"), sep = "+") %>% # Warning: unsorted | this adds the TF footprints counts
      unite("accessible_footprint.width", c("accessible_footprint.width", "footprint.width", "noise_footprint.width"), sep = "+") %>% # Warning: unsorted | this adds the TF footprints counts
      # mutate(accessible_footprint.center = gsub("^NA\\+|\\+NA$", "", accessible_footprint.center), accessible_footprint.width = gsub("^NA\\+|\\+NA$", "", accessible_footprint.width)) %>%
      # mutate(accessible_footprint.center = gsub("\\+NA", "", accessible_footprint.center), accessible_footprint.width = gsub("\\+NA", "", accessible_footprint.width)) %>%
      dplyr::select(cluster.nr, Sample, TFBS.cluster, tot.read.count, read.count, accessible_footprint.center, accessible_footprint.width) %>%
      mutate(center = str_split(accessible_footprint.center, "\\+"), width = str_split(accessible_footprint.width, "\\+")) %>%
      dplyr::select(-accessible_footprint.center, -accessible_footprint.width) %>%
      unnest(c(center, width)) %>%
      filter(center!="NA" & width!="NA") %>%
      mutate(center = as.integer(center), width = as.integer(width)) %>%
      mutate(start = center - width/2, end = center + width/2) %>%
      dplyr::select(-center, -width) %>%
      unite("seqnames", c("Sample", "TFBS.cluster"), sep = ";") %>%
      plyranges::as_granges() %>% 
      sort() %>% 
      plyranges::arrange(cluster.nr) -> acc.granges

  }

  x %>%
    select(Sample, TFBS.cluster, TF, TF.name, footprint.center, footprint.width, cluster.nr.list) %>%
    dplyr::rename(center = "footprint.center", width = "footprint.width") %>%
    mutate(start = center - width/2, end = center + width/2) %>%
    select(-center, -width) %>%
    unite("seqnames", c("Sample", "TFBS.cluster"), sep = ";") %>%
    plyranges::as_granges() -> TFs.granges

  if(!is_null(window.width)){
    TFs.granges %<>% IRanges::resize(., width = window.width, fix = "center")
  }

  TFs.granges %>%
    as_tibble() %>%
    plyranges::mutate(cluster.nr.list = str_split(cluster.nr.list, "\\+")) %>%
    unnest(c(cluster.nr.list)) %>%
    unite("seqnames", "seqnames", "cluster.nr.list", "TF.name", sep = ";") -> tmp

  if(data.type == "F1"){

    plyranges::find_overlaps(TFs.granges, nuc.granges, minoverlap = 2) %>% # this inherently throws away nucleosomes detected on the same partition as TFs
      data.frame() %>%
      select(-c(start, end, width, strand)) %>% # in this line and the one below I use the cluster.nr to make sure I count each partition only once
      distinct() %>%
      unite("seqnames", "seqnames", "cluster.nr", "TF.name", sep = ";") %>%
      filter(!seqnames %in% tmp$seqnames) %>% # this grants that we don't overcount because of partial TF footprints
      separate("seqnames", into = c("Sample", "TFBS.cluster", "cluster.nr", "TF.name"), sep = ";", extra = "merge") %>%
      group_by(Sample, TFBS.cluster, tot.A, tot.R, TF.name) %>%
      summarise(nucleosome.A = sum(A), nucleosome.R = sum(R), cluster.nr.list = paste(unique(cluster.nr), collapse = "+"), .groups = "drop") -> aggregated.nucleosome.counts
    
    aggregated.nucleosome.counts %>%
      mutate(cluster.nr = str_split(cluster.nr.list, "\\+")) %>%
      unnest(c(cluster.nr)) %>%
      unite("seqnames", "Sample", "TFBS.cluster", "cluster.nr", "TF.name", sep = ";") -> tmp.nuc
    
    aggregated.nucleosome.counts %<>%
      dplyr::select(-cluster.nr.list) # to cleanup the return.df
    
    plyranges::find_overlaps(TFs.granges, unr.granges, minoverlap = 2) %>% # this inherently throws away nucleosomes detected on the same partition as TFs
      data.frame() %>%
      select(-c(start, end, width, strand)) %>% # in this line and the one below I use the cluster.nr to make sure I count each partition only once
      distinct() %>%
      unite("seqnames", "seqnames", "cluster.nr", "TF.name", sep = ";") %>%
      filter(!seqnames %in% c(tmp$seqnames, tmp.nuc$seqnames)) %>% # this grants that we don't overcount because of partial TF/nuc footprints
      separate("seqnames", into = c("Sample", "TFBS.cluster", "cluster.nr", "TF.name"), sep = ";", extra = "merge") %>%
      group_by(Sample, TFBS.cluster, tot.A, tot.R, TF.name) %>%
      summarise(unrecognized.A = sum(A), unrecognized.R = sum(R), cluster.nr.list = paste(unique(cluster.nr), collapse = "+"), .groups = "drop") %>%
      mutate(cluster.nr = str_split(cluster.nr.list, "\\+")) %>%
      unnest(c(cluster.nr)) %>%
      unite("seqnames", "Sample", "TFBS.cluster", "cluster.nr", "TF.name", sep = ";") -> tmp.unr
    
    acc.granges %<>% plyranges::group_by(cluster.nr, tot.R, tot.A, R, A) %>% plyranges::reduce_ranges() # Note: this works because TF footprints are included
    Overlaps.acc = IRanges::findOverlaps(TFs.granges, acc.granges, minoverlap = 2, ignore.strand = TRUE) # minoverlap = 2 grants that we don't overcount because of TF-accessibilty boarders
    TFs.granges[queryHits(Overlaps.acc)] %>%
      plyranges::mutate(
        cluster.nr = acc.granges[subjectHits(Overlaps.acc)]$cluster.nr,
        tot.R = acc.granges[subjectHits(Overlaps.acc)]$tot.R,
        tot.A = acc.granges[subjectHits(Overlaps.acc)]$tot.A,
        R = acc.granges[subjectHits(Overlaps.acc)]$R,
        A = acc.granges[subjectHits(Overlaps.acc)]$A,
        acc.width = width(acc.granges[subjectHits(Overlaps.acc)])
      ) %>%
      data.frame() %>%
      select(-c(start, end, width, strand)) %>% # in this line and the one below I use the cluster.nr to make sure I count each partition only once
      distinct() %>%
      unite("seqnames", "seqnames", "cluster.nr", "TF.name", sep = ";") %>%
      filter(!seqnames %in% c(tmp.nuc$seqnames, tmp.unr$seqnames)) %>% # this grants that we don't overcount because of partial nuc/unr footprints. N.b.: I has also tmp$seqnames here
      separate("seqnames", into = c("Sample", "TFBS.cluster", "cluster.nr", "TF.name"), sep = ";", extra = "merge") %>%
      group_by(Sample, TFBS.cluster, tot.A, tot.R, TF.name) %>%
      summarise(accessible.A = sum(A), accessible.R = sum(R), 
                acc.width.distro = list(setNames(acc.width, cluster.nr)), 
                acc.read.count.distro.R = list(setNames(R, cluster.nr)),
                acc.read.count.distro.A = list(setNames(A, cluster.nr)),
                .groups = "drop") -> aggregated.accessible.counts
    
    left_join(left_join(
      utils::type.convert(x, as.is = TRUE),
      aggregated.nucleosome.counts,
      by = c("Sample", "TFBS.cluster", "tot.R", "tot.A", "TF.name")),
      aggregated.accessible.counts,
      by = c("Sample", "TFBS.cluster", "tot.R", "tot.A", "TF.name")
    ) -> return.df
    
    return.df %<>%
      mutate(R = ifelse(is.na(R), 0, R), A = ifelse(is.na(A), 0, A)) %>%
      mutate(nucleosome.R = ifelse(is.na(nucleosome.R), 0, nucleosome.R), nucleosome.A = ifelse(is.na(nucleosome.A), 0, nucleosome.A)) %>%
      mutate(accessible.R = ifelse(is.na(accessible.R), 0, accessible.R), accessible.A = ifelse(is.na(accessible.A), 0, accessible.A))

  } else if (data.type %in% c("PRA", "WT")){

    plyranges::find_overlaps(TFs.granges, nuc.granges, minoverlap = 2) %>% # this inherently throws away nucleosomes detected on the same partition as TFs
      data.frame() %>%
      select(-c(start, end, width, strand)) %>% # in this line and the one below I use the cluster.nr to make sure I count each partition only once
      distinct() %>%
      unite("seqnames", "seqnames", "cluster.nr", "TF.name", sep = ";") %>%
      filter(!seqnames %in% tmp$seqnames) %>% # this grants that we don't overcount because of partial TF footprints
      separate("seqnames", into = c("Sample", "TFBS.cluster", "cluster.nr", "TF.name"), sep = ";", extra = "merge") %>%
      group_by(Sample, TFBS.cluster, tot.read.count, TF.name) %>%
      summarise(nucleosome.read.count = sum(read.count), cluster.nr.list = paste(unique(cluster.nr), collapse = "+"), .groups = "drop") -> aggregated.nucleosome.counts
    
    aggregated.nucleosome.counts %>%
      mutate(cluster.nr = str_split(cluster.nr.list, "\\+")) %>%
      unnest(c(cluster.nr)) %>%
      unite("seqnames", "Sample", "TFBS.cluster", "cluster.nr", "TF.name", sep = ";") -> tmp.nuc
    
    aggregated.nucleosome.counts %<>%
      dplyr::select(-cluster.nr.list) # to cleanup the return.df
    
    plyranges::find_overlaps(TFs.granges, unr.granges, minoverlap = 2) %>% # this inherently throws away nucleosomes detected on the same partition as TFs
      data.frame() %>%
      select(-c(start, end, width, strand)) %>% # in this line and the one below I use the cluster.nr to make sure I count each partition only once
      distinct() %>%
      unite("seqnames", "seqnames", "cluster.nr", "TF.name", sep = ";") %>%
      filter(!seqnames %in% c(tmp$seqnames, tmp.nuc$seqnames)) %>% # this grants that we don't overcount because of partial TF/nuc footprints
      separate("seqnames", into = c("Sample", "TFBS.cluster", "cluster.nr", "TF.name"), sep = ";", extra = "merge") %>%
      group_by(Sample, TFBS.cluster, tot.read.count, TF.name) %>%
      summarise(unrecognized.read.count = sum(read.count), cluster.nr.list = paste(unique(cluster.nr), collapse = "+"), .groups = "drop") %>%
      mutate(cluster.nr = str_split(cluster.nr.list, "\\+")) %>%
      unnest(c(cluster.nr)) %>%
      unite("seqnames", "Sample", "TFBS.cluster", "cluster.nr", "TF.name", sep = ";") -> tmp.unr
    
    acc.granges %<>% plyranges::group_by(cluster.nr, tot.read.count, read.count) %>% plyranges::reduce_ranges() # Note: this works because TF footprints are included
    Overlaps.acc = IRanges::findOverlaps(TFs.granges, acc.granges, minoverlap = 2, ignore.strand = TRUE) # minoverlap = 2 grants that we don't overcount because of TF-accessibilty boarders
    TFs.granges[queryHits(Overlaps.acc)] %>%
      plyranges::mutate(
        cluster.nr = acc.granges[subjectHits(Overlaps.acc)]$cluster.nr,
        tot.read.count = acc.granges[subjectHits(Overlaps.acc)]$tot.read.count,
        read.count = acc.granges[subjectHits(Overlaps.acc)]$read.count,
        acc.width = width(acc.granges[subjectHits(Overlaps.acc)])
      ) %>% 
      data.frame() %>%
      select(-c(start, end, width, strand)) %>% # in this line and the one below I use the cluster.nr to make sure I count each partition only once
      distinct() %>%
      unite("seqnames", "seqnames", "cluster.nr", "TF.name", sep = ";") %>%
      filter(!seqnames %in% c(tmp.nuc$seqnames, tmp.unr$seqnames)) %>% # this grants that we don't overcount because of partial nuc/unr footprints. N.b.: I has also tmp$seqnames here
      separate("seqnames", into = c("Sample", "TFBS.cluster", "cluster.nr", "TF.name"), sep = ";", extra = "merge") %>%
      group_by(Sample, TFBS.cluster, tot.read.count, TF.name) %>%
      summarise(
        accessible.read.count = sum(read.count), 
        acc.width.distro = list(setNames(acc.width, cluster.nr)), 
        acc.read.count.distro = list(setNames(read.count, cluster.nr)), 
        .groups = "drop") -> aggregated.accessible.counts
    
    # aggregated.accessible.counts %>%
    #   unnest(c(acc.width.distro, acc.read.count.distro))
    
    left_join(left_join(
      utils::type.convert(x, as.is = TRUE),
      aggregated.nucleosome.counts,
      by = c("Sample", "TFBS.cluster", "tot.read.count", "TF.name")), 
      aggregated.accessible.counts,
      by = c("Sample", "TFBS.cluster", "tot.read.count", "TF.name")
    ) -> return.df
    
    return.df %<>%
      mutate(read.count = ifelse(is.na(read.count), 0, read.count)) %>%
      mutate(nucleosome.read.count = ifelse(is.na(nucleosome.read.count), 0, nucleosome.read.count)) %>%
      mutate(accessible.read.count = ifelse(is.na(accessible.read.count), 0, accessible.read.count))

  }

  return(return.df)

}

#' this needs ChIP_data_dictionary
add.ChIP.score = function(x, TFBSs, ChIP, genome){
  
  Reduce(rbind, lapply(unique(x$TF), function(tf){
    print(tf)
    if(!is.null(ChIP[[tf]])){
      motifs = TFBSs[filter(x, TF == tf)$TF.name]
      clObj = parallel::makeCluster(16)
      motif.counts = QuasR::qCount(GetQuasRprj(ChIP[[tf]], genome = genome), query = IRanges::resize(motifs, 200, "center"), clObj = clObj)
      parallel::stopCluster(clObj)
      
      left_join(
        filter(x, TF == tf),
        data.frame(TF.name = rownames(motif.counts), ChIP = as.integer(motif.counts[,tf])),
        by = "TF.name"
      ) -> return.df
    } else {
      x %>%
        filter(TF == tf) %>%
        mutate(ChIP = NA) -> return.df
    }
    return(return.df)
  })) -> x_ChIP
  return(x_ChIP)
  
}

#' @param master.table raw
#' @param data.type F1/WT/PRA. WT and PRA have the same effect
#' @param GenomicTiles
#' @param skip.chip T/F whether to skip resource intensive chip-seq quantification
make.chromatin.influence.df = function(master.table, data.type = "F1", GenomicTiles = NULL, skip.chip = FALSE, genome = BSgenome.Mmusculus.UCSC.mm10){
  
  # Not sure where these come from, I don't have time to dig in
  master.table %<>%
    filter(TFBS.cluster != "TRUE")
  
  master.table %>%
    # ANNOTATE CGI
    mutate(CGI = ifelse(TFBS.cluster %in% names(GenomicTiles)[GenomicTiles$CGI], TRUE, FALSE)) %>% 
    # ANNOTATE TILES WITHOUT TF FOOTPRINTS
    annotated_tiles_without_footprints(GenomicTiles = GenomicTiles, TFBSs = TFBSs, data.type = data.type) %>%
    # ADD DELTA PWM
    add.delta.pwm(., TFBSs = TFBSs, data.type = data.type) %>%
    # filter.bound.TFBSs(., TFBSs = TFBSs) %>%
    aggregate.footprints.by.TF.annotation(., TFBSs = TFBSs, data.type = data.type) -> master.table.GenVar
    # filter.TFs.for.binding.frequency(., THR = 0.0, x.type = "singleTF", data.type = data.type) 
  
  # Add nucleosomes
  master.table.GenVar %>%
    mutate(footprint.center = start(IRanges::resize(TFBSs[TF.name], 1, "center")), footprint.width = width(TFBSs[TF.name])) %>% # N.b.: here I restrict the focus on motif CA rather than TF footprint
    add.nucleosome.occupancy(., master.table = master.table, data.type = data.type) -> master.table.GenVar.nuc
  
  # Add Gene Disregulation
  # left_join(
  #   master.table.GenVar.nuc,
  #   TileGeneAssociation %>%
  #     as.data.frame() %>%
  #     select(dist_to_TSS, Status, external_gene_name, Tile) %>%
  #     dplyr::rename(TFBS.cluster = "Tile"),
  #   by = "TFBS.cluster"
  # ) -> master.table.GenVar.nuc.GeneExpr
  master.table.GenVar.nuc -> master.table.GenVar.nuc.GeneExpr
  
  master.table.GenVar.nuc.GeneExpr %>%
    
    mutate(TF = gsub("\\..*$", "", TF)) %>%
    
    filter(!TF %in% c("Nanog", "Nanog-partner", "B-Box")) %>%
    
    mutate(TF_old = TF) %>%
    mutate(
      TF = str_to_title(
        gsub("^MAFK$", "BACH1::MAFK",
             gsub("^SMAD3$", "SMAD2::SMAD3::SMAD4",
                  gsub("POU5F1|POU5F1B", "OCT4", toupper(
                    gsub("_.*$", "", TF_old
                    ))))))) %>%
    # Motif strength
    filter(BL6.absScore >= 10) %>%
    mutate(motif.change = factor(ifelse(delta.pwm/abs(BL6.absScore) > 1.5, "g.o.f.", ifelse(delta.pwm/abs(BL6.absScore) < -.5, "l.o.f.", "no")), levels = c("l.o.f.", "no", "g.o.f."))) %>%
    
    # group_by(TF) %>%
    # {if(data.type == "F1") mutate(., motif.strength.percentile = cut(BL6.absScore, breaks = 10, labels = seq(10, 100, 10))) else mutate(., motif.strength.percentile = NA)} %>%
    # ungroup() %>%
    # {if(data.type == "F1") filter(., (motif.change %in% c("l.o.f.", "no") & motif.strength.percentile %in% seq(80,100,10)) | (motif.change == "g.o.f." & motif.strength.percentile %in% seq(10,50,10))) else .} %>%
    
    # Filter the noise
    {if(data.type == "F1") mutate(., A = ifelse(is.na(A), 0, A), R = ifelse(is.na(R), 0, R)) else mutate(., read.count = ifelse(is.na(read.count), 0, read.count))} %>% 
    # {if(data.type == "F1") mutate(., accessible.A = accessible.A + A, accessible.R = accessible.R + R) else mutate(., accessible.read.count = accessible.read.count + read.count)} %>% ####### Now the following line happens inside add.nucleosome.occupancy
    
    # Prior on accessibility
    # {if(data.type == "F1") filter(., accessible.R/tot.R >= 0.4) else if(data.type == "WT") filter(., accessible.read.count/tot.read.count >= 0.4) else if(data.type == "PRA") .} %>%
    
    {if(data.type == "F1") mutate(., nucleosome.delta = (nucleosome.A/tot.A)-(nucleosome.R/tot.R)) else mutate(., nucleosome.delta = NA)} %>%
    {if(data.type == "F1") mutate(., accessible.delta = (accessible.A/tot.A)-(accessible.R/tot.R)) else mutate(., accessible.delta = NA)} %>%
    gather(feature, delta, nucleosome.delta, accessible.delta) %>%
    mutate(feature = gsub("\\.delta", "", feature)) %>%
    
    filter(feature == "accessible") %>%
    group_by(TF) %>%
    {if(data.type == "F1" & nrow(.)>1) filter(., sum(motif.change == "l.o.f.") >= 1) else {.}} %>% # N.b.: sum(motif.change == "l.o.f.") >= 1 is to filter out TFs with no perturbations when doing GW analysis. nrow(.)>1 is to let the dataframe pass when plotting single sites.
    ungroup() -> chromatin.influence.df
  
  # REDUCE REDUNDANCY
  #        1. GenomicTiles redundancy: pick the tile which centers the TFBS the most
  chromatin.influence.df %<>% 
    mutate(
      tile_center_coord.tmp = start(IRanges::resize(GenomicTiles[TFBS.cluster], width = 1, fix = "center")), 
      TFBS_center_coord.tmp = start(IRanges::resize(TFBSs[TF.name], width = 1, fix = "center")),
      tile_TFBS_centers_dist.tmp = abs(tile_center_coord.tmp - TFBS_center_coord.tmp)
    ) %>%
    group_by(Sample, TF.name) %>%
    slice_min(tile_TFBS_centers_dist.tmp) %>%
    ungroup() %>%
    select(-c(tile_center_coord.tmp, TFBS_center_coord.tmp, tile_TFBS_centers_dist.tmp))
  #        2. F1s redundancy         : when the two F1s present the same motif change, pick the one with the best SMF coverage
  if(data.type == "F1"){
    chromatin.influence.df %<>%
      arrange(TF.name) %>%
      group_by(TF.name, delta.pwm) %>%
      slice_max(tot.A+tot.R) %>%
      ungroup()
  }
  #        3. TFBSs redundancy       : when multiple PWMs are used for the same TF, pick the one with the most extreme delta pwm
  TFBSs.tmp = TFBSs[chromatin.influence.df$TF.name]
  TFBSs.tmp$reduced.TF = chromatin.influence.df$TF
  GenomicRanges::findOverlaps(TFBSs.tmp, TFBSs.tmp, minoverlap = 5, ignore.strand = TRUE) -> ov
  g = igraph::graph_from_adj_list(as.list(ov), mode = "all")
  comp = igraph::components(g)
  sets = igraph::groups(comp)
  names(sets) = paste0("idx.", names(sets))
  data.frame(
    grouping.factor = rep(names(sets), lengths(sets)), idx = unlist(sets, use.names = FALSE)
  ) %>% arrange(idx) -> grouping.factor.tmp
  chromatin.influence.df %<>%
    mutate(grouping.factor = grouping.factor.tmp$grouping.factor) %>%
    group_by(Sample, TF, grouping.factor) %>%
    mutate(hemi = ifelse(str_detect(TF_old, "_h$|_H$"), TRUE, FALSE)) %>%
    arrange(BL6.absScore) %>%
    slice_head(n=1) %>%
    ungroup() %>%
    select(-c(grouping.factor, hemi))
  
  chromatin.influence.df %<>%
    mutate(fishers.l2OR = NA, fishers.pval = NA)
    
  # Add ChIP
  if(isFALSE(skip.chip)){
    chromatin.influence.df %<>% add.ChIP.score(., TFBSs = TFBSs, ChIP = ChIP_data_dictionary, genome = genome)
  }
  
  return(chromatin.influence.df)
  
}

#' adds CRE metadata & reduces metrics
#' @param x output of make.chromatin.influence.df
#' @param cre.annotation chromHMM
#' @param chip.thr ChIP_thresholds_dictionary_lenient
process.CA.df = function(x, cre.annotation, chip.thr){
  
  x %>%
    left_join(., cre.annotation, by = "TFBS.cluster") %>%
    left_join(., chip.thr, by = "TF") %>%
    mutate(ChIP_annotation = factor(case_when(
      ChIP == 0 ~ "unbound",
      ChIP <= ChIP_threshold ~ "weak",
      ChIP > ChIP_threshold ~ "bound"
    ), levels = c("unbound", "weak", "bound"))) %>%
    mutate(width = unlist(pmap(list(acc.width.distro, acc.read.count.distro), function(x, w){weighted.mean(x, w)}))) %>%
    mutate(CA = (accessible.read.count/tot.read.count)*100) %>%
    unnest(c(acc.read.count.distro, acc.width.distro)) %>%
    filter(acc.width.distro >= 100) %>%
    group_by(Sample, TFBS.cluster, TF.name, TF, chromHMM_annotation, ChIP_annotation, tot.read.count, width, CA) %>%
    summarise(
      CA_regulatory_count = sum(acc.read.count.distro),
      acc.width.distro = list(acc.width.distro), 
      acc.read.count.distro = list(acc.read.count.distro), 
      .groups = "drop") %>%
    mutate(CA_regulatory = CA_regulatory_count / tot.read.count * 100) %>%
    mutate(
      width_regulatory = unlist(pmap(list(acc.width.distro, acc.read.count.distro), function(x, w){weighted.mean(x, w)})),
    ) %>%
    dplyr::select(-c(acc.width.distro, acc.read.count.distro)) -> df
  
  return(df)
  
}

#' adds CRE metadata & reduces metrics
#' @param x output of make.chromatin.influence.df
#' @param cre.annotation chromHMM
#' @param chip.thr ChIP_thresholds_dictionary_lenient
process.CA.df_f1 = function(x, cre.annotation, chip.thr){
  
  cluster <- multidplyr::new_cluster(16)
  x %>%
    filter(motif.change != "g.o.f.") %>%
    left_join(., cre.annotation, by = "TFBS.cluster") %>%
    left_join(., chip.thr, by = "TF") %>%
    mutate(ChIP_annotation = factor(case_when(
      ChIP == 0 ~ "unbound",
      ChIP <= ChIP_threshold ~ "weak",
      ChIP > ChIP_threshold ~ "bound"
    ), levels = c("unbound", "weak", "bound"))) %>%
    mutate(
      width_R = unlist(pmap(list(acc.width.distro, acc.read.count.distro.R), function(x, w){weighted.mean(x, w)})),
      width_A = unlist(pmap(list(acc.width.distro, acc.read.count.distro.A), function(x, w){weighted.mean(x, w)}))
      ) %>%
    mutate(CA_R = (accessible.R/tot.R)*100, CA_A = (accessible.A/tot.A)*100) %>%
    unnest(c(acc.read.count.distro.R, acc.read.count.distro.A, acc.width.distro)) %>%
    filter(acc.width.distro >= 100) %>%
    group_by(Sample, TFBS.cluster, TF.name, TF, motif.change, chromHMM_annotation, ChIP_annotation, tot.R, tot.A, width_R, width_A, CA_R, CA_A) %>%
    summarise(
      CA_regulatory_count_R = sum(acc.read.count.distro.R), CA_regulatory_count_A = sum(acc.read.count.distro.A), 
      acc.width.distro = list(acc.width.distro), 
      acc.read.count.distro.R = list(acc.read.count.distro.R), 
      acc.read.count.distro.A = list(acc.read.count.distro.A),
      .groups = "drop") %>%
    mutate(CA_regulatory_R = CA_regulatory_count_R / tot.R * 100, CA_regulatory_A = CA_regulatory_count_A / tot.A * 100) %>%
    mutate(
      width_regulatory_R = unlist(pmap(list(acc.width.distro, acc.read.count.distro.R), function(x, w){weighted.mean(x, w)})),
      width_regulatory_A = unlist(pmap(list(acc.width.distro, acc.read.count.distro.A), function(x, w){weighted.mean(x, w)}))
    ) %>%
    dplyr::select(-c(acc.width.distro, acc.read.count.distro.R, acc.read.count.distro.A)) %>%
    rowwise() %>%
    multidplyr::partition(cluster) %>%
    mutate(pval = fisher.test(matrix(c(CA_regulatory_count_A,tot.A-CA_regulatory_count_A,CA_regulatory_count_R,tot.R-CA_regulatory_count_R),ncol = 2), alternative = "less")$p.value) %>%
    ungroup() %>%
    collect() %>%
    mutate(CA.delta = CA_regulatory_A - CA_regulatory_R) -> df
  remove(cluster)
  df$p.adj = NA
  df$p.adj[df$motif.change=="l.o.f." & df$ChIP_annotation=="bound" & !is.na(df$ChIP_annotation)] = p.adjust(p = filter(df, motif.change=="l.o.f." & ChIP_annotation=="bound")$pval, method = "BH")
  df %<>% mutate(change = ifelse(p.adj <= 0.05 & CA.delta < 0, "down", ifelse(p.adj <= 0.05 & CA.delta > 0, "up", "no")))
  
  return(df)
  
}

