#' This is useful for bait capture or genome-wide data when
#' wanting to cluster all covered sites before knowing which ones are covered
#' 
#' @param sampleSheet
#' @param Samples
#' @param pool.samples Whether to pool or not samples for clustering & footprints detection. Defaults to TRUE. 
#' @param genome
#' @param TFBS.clusters object resulting from a Arrange_TFBSs_clusters call
#' @param SNPs all SNPs (both species)
#' @param TFBSs Granges with all the annotated TFBSs
#' @param CytosinesToMask
#' @param species.string.match Used to pick SNPs for plotting. Should be a match for how the Cast species is encoded in the SampleName string. Everything that doesn't match will be considered Spretus
#' @param deduplicate TRUE or FALSE (default). Useful with samples with high PCR duplication
#' @param params bait.capture.parameters.list
#' @param cores
#' 
Unsupervised.clustering.MultiSite.wrapper = function(
  sampleSheet, Samples, pool.samples = TRUE, genome, TFBS.clusters, # passed to - TFBS.clusters.coverage -
  SNPs, MethylationCallingWindows, TFBSs, CytosinesToMask, species.string.match = "", 
  deduplicate = FALSE,
  params = NULL, cores = 40, tryCatchLog.write.error.dump.file = FALSE
  ){
  
  # MethylationCallingWindows = SingleMoleculeFootprinting::Create_MethylationCallingWindows(
  #   TFBS_cluster_coordinates = TFBS.clusters$ClusterCoordinates, 
  #   genomic.seqlenghts = GenomeInfoDb::seqlengths(genome), 
  #   fix.window.size = TRUE, max.window.size = 50)
 
  print(paste0("Analysing ", length(MethylationCallingWindows), " methylation calling windows"))
  
  parallel::mclapply(seq_along(MethylationCallingWindows), function(i){
    
    print(paste0("STARTING ANALYSIS OF METHYLATION CALLING WINDOW ", i, "/", length(MethylationCallingWindows)))
    
    # Unpack parameters ####
    CallContextMethylation.coverage = params$CallContextMethylation[["coverage"]]
    CallContextMethylation.ConvRate.thr = params$CallContextMethylation[["ConvRate.thr"]]
    ########################
    
    # tryCatchLog::tryCatchLog({
    tryCatch({
      
      CurrentWindow = MethylationCallingWindows[i]
      CurrentWindow_ext = CurrentWindow
      start(CurrentWindow_ext) = start(CurrentWindow_ext) - 250
      end(CurrentWindow_ext) = end(CurrentWindow_ext) + 250
      if(start(CurrentWindow_ext) < 1){start(CurrentWindow_ext) = 1}
      if(end(CurrentWindow_ext) > seqlengths(genome)[[as.character(seqnames(CurrentWindow_ext))]]){end(CurrentWindow_ext) = seqlengths(genome)[[as.character(seqnames(CurrentWindow_ext))]]}
      
      Experiment = SingleMoleculeFootprinting::DetectExperimentType(Samples = Samples)
      CallContextMethylation(sampleSheet = sampleSheet, sample = Samples,
                             genome = genome, RegionOfInterest = CurrentWindow_ext,
                             coverage = CallContextMethylation.coverage, ConvRate.thr = CallContextMethylation.ConvRate.thr,
                             returnSM = TRUE, clObj = NULL) -> Methylation
      
      if (!is.null(CytosinesToMask)){
        if(pool.samples){
          Methylation <- MaskSNPs(Methylation = Methylation, CytosinesToMask = CytosinesToMask,
                                  MaskSMmat = TRUE, Experiment = Experiment)
        } else {
          Methylation <- MaskSNPs2(Methylation = Methylation, CytosinesToMaks = CytosinesToMask,
                                   MaskSMmat = TRUE, Experiment = Experiment)
        }
      }
      
      # NOTE: for now let's just discard endogenous CpG info, as I don't see myself using it in this prj
      if (Experiment == "NO"){
        Methylation[[1]] = Methylation[[1]]$DGCHN
        Methylation[[2]] = lapply(Methylation[[2]], function(x){x$DGCHN})
      }
      
      if (length(Methylation[[1]]) == 0){
        stop("No cytosines found covering the region")
      }
      
      TFBS.clusters.in.MethylationCallingWindows = list()
      TFBS.clusters.in.MethylationCallingWindows$ClusterCoordinates = plyranges::filter_by_overlaps(TFBS.clusters$ClusterCoordinates, CurrentWindow, minoverlap = min(width(TFBS.clusters$ClusterCoordinates)))
      
      if (pool.samples){
        message("Pooling samples before clustering...iterating over TFBS.clusters")
        iteration.list = lapply(names(TFBS.clusters.in.MethylationCallingWindows$ClusterCoordinates), function(i){Samples})
        names(iteration.list) = names(TFBS.clusters.in.MethylationCallingWindows$ClusterCoordinates)
      } else {
        message("Not pooling samples before clustering...iterating over TFBS.clusters & samples")
        iteration.list = split(rep(Samples, length(TFBS.clusters.in.MethylationCallingWindows$ClusterCoordinates)), seq(length(Samples)*length(TFBS.clusters.in.MethylationCallingWindows$ClusterCoordinates)))
        names(iteration.list) = rep(names(TFBS.clusters.in.MethylationCallingWindows$ClusterCoordinates), each = length(Samples))
      }
      
      lapply(seq_along(iteration.list), function(j){

        # Unpack parameters ####
        output.params = params$output.params
        ########################
        
        cluster = names(iteration.list)[j]
        current.sample = iteration.list[[j]]
        RegionOfInterest = TFBS.clusters.in.MethylationCallingWindows$ClusterCoordinates[cluster]
        RegionOfInterest_ext = resize(RegionOfInterest, width = 500, fix = "center")
        TFBSs.overlapping.locus = plyranges::filter_by_overlaps(TFBSs[seqnames(TFBSs) == seqnames(RegionOfInterest_ext)], RegionOfInterest_ext)

        print(paste0("iteration list step ", j,"/", length(iteration.list), " [MethylationCallingWindow ", i, "/", length(MethylationCallingWindows), "]: ", names(RegionOfInterest)))
        message(paste0("iteration list step ", j,"/", length(iteration.list), " [MethylationCallingWindow ", i, "/", length(MethylationCallingWindows), "]: ", names(RegionOfInterest)))

        #tryCatchLog::tryCatchLog({
        tryCatch({
          
          read.origin = lapply(Methylation[[2]], rownames)
          binary.matrix = Reduce(SingleMoleculeFootprinting::rbind.fill.matrix.sparse, Methylation[[2]])
          binary.matrix = binary.matrix[,colnames(binary.matrix) %in% as.character(start(RegionOfInterest_ext):end(RegionOfInterest_ext)), drop=FALSE]
          
          if(deduplicate){
            message("removing PCR duplicates and incomplete molecules")
            incomplete.molecules.idx = as.logical(rowSums(binary.matrix == 0) > 0)
            binary.matrix = binary.matrix[!incomplete.molecules.idx,]
            duplicated.molecules.idx = as.logical(duplicated(as.matrix(binary.matrix)))
            binary.matrix = binary.matrix[!duplicated.molecules.idx,]
          }
          
          if (ncol(binary.matrix) > 0){

            compute.pam.clustering(
              binary.matrix = binary.matrix, TFBSs.overlapping.locus = TFBSs.overlapping.locus, 
              RegionOfInterest = RegionOfInterest, RegionOfInterest_ext = RegionOfInterest_ext,
              params = params, read.origin = read.origin, output.params = output.params
              ) -> res
            
            message(paste0("> pam completed: ", names(RegionOfInterest)))
            
            Reduce(rbind,
                   lapply(seq_along(read.origin), function(z){
                     
                     current.sample = names(read.origin)[z]
                     current.reads = read.origin[[z]]
                     all.cluster.nrs = seq(params$pam[["k"]])
                     clustering.reads = res$pam.res[names(res$pam.res) %in% current.reads]
                     
                     if(length(clustering.reads) > 0){
                       
                       clustering.reads %>%
                         table() %>%
                         as.data.frame() %>%
                         dplyr::rename("cluster.nr" = ".", "read.count" = "Freq") %>%
                         mutate(Sample = current.sample, TFBS.cluster = names(RegionOfInterest)) -> df
                       
                       # expand unused levels
                       if (length(levels(df$cluster.nr)) != length(all.cluster.nrs)){
                         missing.cluster.nr = all.cluster.nrs[!all.cluster.nrs %in% df$cluster.nr]
                         df = rbind(
                           df, 
                           data.frame(cluster.nr = as.factor(missing.cluster.nr), read.count = 0, Sample = current.sample, TFBS.cluster = names(RegionOfInterest))
                         )
                       }
                       
                     } else {
                       
                       data.frame(matrix(nrow = 0, ncol = 4, dimnames = list(NULL, c("cluster.nr", "read.count", "Sample", "TFBS.cluster")))) -> df
                       
                     }
                     
                     return(df)
                     
                   })) -> res$cluster.coverage.df
            message("> collected coverage.df")
            res$read.origin = NULL
            res$RegionOfInterest = RegionOfInterest
            res$footprints.df %<>%
              mutate(
                TFBS.cluster = names(RegionOfInterest)#,
              )
            
            return(res)
            
          } else {
            NULL
          }
          
        },
        error = function(e){NULL}#,
        # write.error.dump.file = tryCatchLog.write.error.dump.file, 
        # write.error.dump.folder = "/scratch/barzaghi/tryCatchLog_dump/"
        ) -> res
        
        if (!is.null(res)){
          
          return(list(
            # 1.
            # distance.df = NULL,
            # 2.
            footprints.df = res$footprints.df,
            # 3.
            cluster.coverage.df = res$cluster.coverage.df
          ))
          
        } else {
          return(list(
            # distance.df = NULL, 
            footprints.df = NULL, 
            cluster.coverage.df = NULL))
        }
        
      }) -> all.iteration.list.res
      
      print(paste0("Succesfully analysed ", length(all.iteration.list.res), "/", length(iteration.list), " in MethylationCallingWindow ", i, "/", length(MethylationCallingWindows)))
      
      return(
        list(
          # unsupervised.distance.results.df = NULL,
          unsupervised.footprints.detection.results.df = data.table::rbindlist(map(all.iteration.list.res, "footprints.df")),
          cluster.coverage.results.df = data.table::rbindlist(map(all.iteration.list.res, "cluster.coverage.df"))
        )
      )
      
    }, 
    error = function(e){NULL}#,
    # write.error.dump.file = tryCatchLog.write.error.dump.file, 
    # write.error.dump.folder = "/scratch/barzaghi/tryCatchLog_dump/"
    ) -> intermediate.return.object
    
    return(intermediate.return.object)
    
    }, mc.cores = ifelse(cores > length(MethylationCallingWindows), length(MethylationCallingWindows), cores), mc.preschedule = FALSE) -> All.MethylationCallingWindows.results
    
  print(paste0("Succesfully analysed ", length(All.MethylationCallingWindows.results), "/", length(MethylationCallingWindows), " MethylationCallingWindows"))
  
  return(
    list(
        # arguments = NULL,
        # info = NULL,
        # unsupervised.distance.results.df = NULL,
        unsupervised.footprints.detection.results.df = data.table::rbindlist(map(All.MethylationCallingWindows.results, "unsupervised.footprints.detection.results.df")), #
        cluster.coverage.results.df = data.table::rbindlist(map(All.MethylationCallingWindows.results, "cluster.coverage.results.df")) 
        )
  )
  
}














