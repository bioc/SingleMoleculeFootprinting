#' Consutruct column annotation df for pam.results.heatmap
#'
# Heatmap.col.annotation.df = function(sliding.window.matrix, RegionOfInterest, TFBSs = NULL, PlottingSNPs = NULL){
#   
#   current.TFBSs = TFBSs[queryHits(findOverlaps(TFBSs, RegionOfInterest))]
#   TBFSs.annotation = rep(NA, ncol(sliding.window.matrix))
#   findOverlaps(current.TFBSs, GRanges(seqnames = unique(seqnames(current.TFBSs)), IRanges(sliding.window.matrix %>% colnames() %>% as.integer(), width = 1))) -> ov
#   TBFSs.annotation[subjectHits(ov)] = as.character(current.TFBSs$TF[queryHits(ov)])
#   TBFSs.annotation.text = rep(NA, ncol(sliding.window.matrix))
#   findOverlaps(IRanges::resize(current.TFBSs, 1, "center"), GRanges(seqnames = unique(seqnames(current.TFBSs)), IRanges(sliding.window.matrix %>% colnames() %>% as.integer(), width = 1))) -> ov
#   TBFSs.annotation.text[subjectHits(ov)] = as.character(current.TFBSs$TF[queryHits(ov)])
#   annotation.df = data.frame(TFBSs = TBFSs.annotation, TFBSs.text = TBFSs.annotation.text)
#   rownames(annotation.df) = colnames(sliding.window.matrix)
#   if(!is.null(PlottingSNPs)){
#     # PlottingSNPs[queryHits(findOverlaps(PlottingSNPs, RegionOfInterest))] %>%
#     #   as_tibble() %>%
#     #   select(start, R, A) %>%
#     #   mutate(start = as.character(start)) %>%
#     #   dplyr::left_join(rownames_to_column(annotation.df, "start"), ., by="start") %>%
#     #   replace(is.na(.), "") %>%
#     #   column_to_rownames("start") -> annotation.df
#     PlottingSNPs[queryHits(findOverlaps(PlottingSNPs, RegionOfInterest))] %>%
#       as_tibble() %>%
#       dplyr::select(start) -> tmp
#     if(nrow(tmp) > 0){
#       data.frame(start = as.character(sort(unlist(c(tmp$start, tmp$start + 1, tmp$start - 1)))), SNPs = "TRUE") %>%
#         distinct() %>%
#         dplyr::left_join(rownames_to_column(annotation.df, "start"), ., by="start") %>%
#         replace(is.na(.), "") %>%
#         column_to_rownames("start") -> annotation.df
#     }
#   }
#   return(annotation.df)
#   
# }

# viz.fix.binary.heatmap = function(mat, final.size = 7){
#   
#   data.columns = which(!colSums(is.na((mat))) == nrow(mat))
#   for(i in seq_along(data.columns)){
#     columns.to.fill = (data.columns[i]-(floor(final.size/2))):(data.columns[i]+(floor(final.size/2)))
#     if(min(columns.to.fill) < 1){
#       columns.to.fill = seq(max(columns.to.fill))
#     }
#     if(max(columns.to.fill) > ncol(mat)){
#       columns.to.fill = min(columns.to.fill):ncol(mat)
#     }
#     mat[,(data.columns[i]-(floor(final.size/2))):(data.columns[i]+(floor(final.size/2)))] =
#       matrix(rep(mat[,data.columns[i],drop=FALSE], final.size), ncol = final.size)
#   }
#   return(mat)
#   
# }

#' #' @import ComplexHeatmap
#' library(ComplexHeatmap)
#' ht_opt$message = FALSE
#' pam.results.heatmap = function(sliding.window.matrix, binary.mat, RegionOfInterest, TFBSs, params, clustering.results, silhouette.results = NULL, footprints.df, main = ""){
#'   
#'   # Unpack parameters ####
#'   PlottingSNPs = params$PlotSingleSiteSMF[["PlottingSNPs"]]
#'   ########################
#'   
#'   if (!is.null(silhouette.results)){
#'     clusters.viz.order = arrange.clusters.by.footprints(footprints.df = footprints.df)
#'     silhouette.results[,] %>% as_tibble() %>% 
#'       mutate(indexes = 1:nrow(.), cluster = factor(cluster, levels = clusters.viz.order)) %>% 
#'       arrange(cluster, desc(sil_width)) -> sorted.silhouette.results
#'     clustering.results = factor(clustering.results[sorted.silhouette.results$indexes], levels = clusters.viz.order)
#'     pl.mat = 1 - sliding.window.matrix[sorted.silhouette.results$indexes,,drop=FALSE]
#'     pl.mat.raw = 1 - make.full.coordinate.matrix(make.SM.mat.dense(binary.mat), params = params)[sorted.silhouette.results$indexes,,drop=FALSE]
#'     left.annotation = HeatmapAnnotation(silhouette = anno_barplot(x = sorted.silhouette.results$sil_width, 
#'                                                                   baseline = 0, ylim = c(-0.5,1),
#'                                                                   gp = gpar(fill = clustering.results, 
#'                                                                             col = clustering.results)), 
#'                                         which = "row", 
#'                                         annotation_name_rot = 90) 
#'   } else {
#'     sorting.index = order(clustering.results)
#'     clustering.results = clustering.results[sorting.index]
#'     pl.mat = 1 - sliding.window.matrix[sorting.index,,drop=FALSE]
#'     
#'     left.annotation = HeatmapAnnotation(cluster = clustering.results, 
#'                                         which = "row", 
#'                                         annotation_name_rot = 90) 
#'   }
#'   
#'   # x.axis label
#'   labels.cols = rep("", ncol(pl.mat))
#'   labels.cols[ncol(pl.mat) %>% seq %>% quantile(., seq(0,1,0.1)) %>% ceiling %>% as.integer] = 
#'     colnames(pl.mat)[ncol(pl.mat) %>% seq %>% quantile(.,seq(0,1,0.1)) %>% ceiling %>% as.integer]
#'   
#'   # colors
#'   if (length(unique(pl.mat)) == 2){
#'     col_func_SMF = structure(2:1, names = 0:1)
#'   } else {
#'     col_func_SMF = circlize::colorRamp2(breaks = seq(0,1,0.01), 
#'                                         colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdYlBu")))(101))
#'   }
#'   
#'   # top annotation
#'   if(is.null(TFBSs)){
#'     top.annotation = NULL
#'   } else if (is.null(PlottingSNPs) & !is.null(TFBSs)) {
#'     col.annotation.df = Heatmap.col.annotation.df(sliding.window.matrix = sliding.window.matrix, RegionOfInterest = RegionOfInterest, TFBSs = TFBSs, PlottingSNPs = PlottingSNPs)
#'     col.annotation.df$TFBSs.text[is.na(col.annotation.df$TFBSs.text)] = ""
#'     TFBSs.colors = stats::setNames(rep("dimgray", length(purrr::discard(unique(col.annotation.df$TFBSs), is.na))), purrr::discard(unique(col.annotation.df$TFBSs), is.na))
#'     if(is.null(PlottingSNPs)){
#'       top.annotation = columnAnnotation(TFBSs.text = anno_text(col.annotation.df$TFBSs.text, rot = 0, just = "centre", gp = gpar(fontsize = 8)),
#'                                         TFBSs = col.annotation.df$TFBSs,
#'                                         col = list(TFBSs = TFBSs.colors),
#'                                         na_col = "white")
#'   }} else if (!is.null(TFBSs) & !is.null(PlottingSNPs)){
#'     col.annotation.df = Heatmap.col.annotation.df(sliding.window.matrix = sliding.window.matrix, RegionOfInterest = RegionOfInterest, TFBSs = TFBSs, PlottingSNPs = PlottingSNPs)
#'     col.annotation.df$TFBSs.text[is.na(col.annotation.df$TFBSs.text)] = ""
#'     TFBSs.colors = stats::setNames(rep("dimgray", length(purrr::discard(unique(col.annotation.df$TFBSs), is.na))), purrr::discard(unique(col.annotation.df$TFBSs), is.na))
#'     SNPs.colors = setNames(rep("dimgray", length(unique(col.annotation.df$SNPs))), unique(col.annotation.df$SNPs))
#'     top.annotation = columnAnnotation(TFBSs.text = anno_text(col.annotation.df$TFBSs.text, rot = 0, just = "centre", gp = gpar(fontsize = 8)),
#'                                       TFBSs = col.annotation.df$TFBSs,
#'                                       # Ref = anno_text(col.annotation.df$R, rot = 0),
#'                                       # Alt = anno_text(col.annotation.df$A, rot = 0),
#'                                       SNPs = col.annotation.df$SNPs,
#'                                       col = list(TFBSs = TFBSs.colors, SNPs = SNPs.colors),
#'                                       na_col = "white")
#'   } else if (is.null(TFBSs) & !is.null(PlottingSNPs)){
#'     stop("I haven't coded this case")
#'   }
#'   
#'   # right annotation
#'   # clusters.table = as.integer(table(sorted.silhouette.results$cluster))
#'   # right.annotation.values = list(Tf.fp = rep(footprints.df$Nr.TF.fp, clusters.table), Nucleosome.fp = rep(footprints.df$Nr.Nucleosome.fp, clusters.table))
#'   # right.annotation = rowAnnotation(footprins = anno_barplot(x = right.annotation.values, gp = gpar(fill = 2:3, col = 2:3), 
#'   #                                                           add_numbers = TRUE, height = unit(1, "cm")), annotation_name_rot = 90)
#'   
#'   left.annotation +
#'     Heatmap(matrix = as.matrix(pl.mat), name = "SMF.sw", column_title = main, column_title_gp = gpar(fontsize = 8),
#'             col = col_func_SMF, width = 4,
#'             cluster_rows = FALSE, cluster_columns = FALSE, 
#'             show_column_names = TRUE, show_row_names = FALSE, 
#'             column_labels = labels.cols, 
#'             row_split = clustering.results, 
#'             border = TRUE, row_gap = unit(1, "mm"),
#'             top_annotation = top.annotation) +
#'     Heatmap(matrix = viz.fix.binary.heatmap(as.matrix(pl.mat.raw), final.size = 7), 
#'             name = "SMF.raw", column_title = "_", 
#'             col = structure(2:1, names = 0:1), width = 4,
#'             cluster_rows = FALSE, cluster_columns = FALSE,
#'             show_column_names = FALSE, show_row_names = FALSE,
#'             column_labels = labels.cols
#'             )# +
#'     # right.annotation
#'   
#' }

#' @param TFBSs.overlapping.locus is assumed to be already subset for the TFBSs of interest
#' @param params list of parameters to pass to lower functions
#' @param output.params list of parameters determining what is returned by the function
#' @param read.origin should be either NULL or a list of 2 vectors containing read names oroginating from the 2 different alleles
#' 
#' @importFrom parallelDist parDist
#' @importFrom cluster pam silhouette
#' @importFrom SingleMoleculeFootprinting PlotSM
compute.pam.clustering = function(binary.matrix, TFBSs.overlapping.locus, RegionOfInterest, RegionOfInterest_ext, 
                                  params, 
                                  read.origin=NULL,
                                  # scan.orphan.footprints,
                                  output.params = output.params){
  
  # Unpack parameters ####
  plot.SW.heatmap = output.params$plot.SW.heatmap
  plot.binary.heatmap = output.params$plot.binary.heatmap
  compute.silhouette = output.params$compute.silhouette
  return.data = output.params$return.data
  
  unsupervised.clustering.coverage = params$unsupervised.clustering.coverage
  max.nr.reads = params$compute.pam.clustering[["max.nr.reads"]]
  method = params$parDist[["method"]]
  threads = as.integer(params$parDist[["threads"]])
  k = params$pam[["k"]]
  ########################
  
  if (!is.null(max.nr.reads)){
    if(nrow(binary.matrix) > max.nr.reads){
      message(paste0("0. Subsetting SW.mat to ", max.nr.reads, " reads"))
      binary.matrix = binary.matrix[sort(sample(seq(nrow(binary.matrix)), max.nr.reads, replace = FALSE)),,drop=FALSE]
    }
  }
  
  message("1. Computing sliding windows")
  SW.mat = compute.sliding.windows(binary.matrix = binary.matrix, TFBSs.to.cover = RegionOfInterest, params = params)
  # try not to smooth
  # SW.mat = binary.matrix %>% make.SM.mat.dense() %>% filter.single.molecule.matrix(., ranges.to.cover = RegionOfInterest, params = params)
  
  if(nrow(SW.mat) < unsupervised.clustering.coverage){
    stop(paste0("The site is covered by less than ", unsupervised.clustering.coverage, " continuous reads...quitting"))
  }
  
  message("2. computing distance matrix")
  distance.matrix = parallelDist::parDist(x = SW.mat, method = method, threads = threads)
  
  message("3. Clustering")
  pam.res = cluster::pam(x = distance.matrix, k = k, diss = TRUE, cluster.only = TRUE, keep.diss = FALSE, keep.data = FALSE)
  new.k = k
  while(any(table(pam.res) < params$detect.footprint[["cytosine.min.coverage"]])){
    message(paste0("partitions too slim detected...retrying with k=", new.k-1))
    new.k = new.k - 1
    if(new.k == 1){stop('Cannot be clustered given the params$detect.footprint[["cytosine.min.coverage"]] required')}
    pam.res = cluster::pam(x = distance.matrix, k = new.k, diss = TRUE, cluster.only = TRUE, keep.diss = FALSE, keep.data = FALSE)
  }
  
  if(compute.silhouette){
    sil.res = cluster::silhouette(x = pam.res, dist = distance.matrix)
  } else {sil.res = NULL}
  
  remove(distance.matrix)
  
  message("4. Recomputing SW mat over full reads")
  # SW.mat = compute.full.read.sliding.windows(binary.matrix = binary.matrix, params = params)
  # SW.mat = SW.mat[names(pam.res),]
  
  message("5. Footprints detection & annotation")
  detect.footprint.wrapper(
    single.molecule.matrix = binary.matrix, clustering.results = pam.res, params = params
    ) -> footprints.df
  footprint.detection.diagnostic.plot(
    single.molecule.matrix = binary.matrix, clustering.results = pam.res, footprints.df = footprints.df, assigned.TFBSs = TFBSs.overlapping.locus
    ) -> footprint.plot
  
  Annotate.footprints.wrapper(
    footprints.df = footprints.df, chromosome = as.character(seqnames(RegionOfInterest)), TFBSs.overlapping.locus = TFBSs.overlapping.locus, 
    RNAseq = NULL, 
    params = params) -> footprints.annotation.results
  
  TFBSs.overlapping.locus = footprints.annotation.results$TFBSs.overlapping.locus
  footprints.df = gather.equivalent.footprints(footprints.df = footprints.annotation.results$footprints.df)
 
  message("6. plotting")
  # if (plot.SW.heatmap){
  #   if(is.null(read.origin)){ # I think this condition is never TRUE, it works fine though
  #     # pam.results.heatmap(sliding.window.matrix = SW.mat, RegionOfInterest = RegionOfInterest, TFBSs = TFBSs.overlapping.locus, params = params,
  #     #                     clustering.results = pam.res, silhouette.results = sil.res, footprints.df = footprints.df,
  #     #                     main = paste(method, "| k=", k, "| ", nrow(sil.res[,]), " reads")) -> pam.heatmap
  #     # results$pam.heatmaps = pam.heatmap
  #   } else {
  #     lapply(seq_along(read.origin), function(i){
  #       reads.to.plot = rownames(SW.mat)[rownames(SW.mat) %in% read.origin[[i]]]
  #       if (length(reads.to.plot) > 0){
  #         reads.to.plot.index = which(rownames(SW.mat) %in% read.origin[[i]])
  #         binary.mat = binary.matrix[reads.to.plot,]
  #         pam.results.heatmap(sliding.window.matrix = SW.mat[reads.to.plot,], binary.mat = binary.mat, RegionOfInterest = RegionOfInterest, TFBSs = TFBSs.overlapping.locus, params = params,
  #                             clustering.results = pam.res[reads.to.plot], silhouette.results = sil.res[reads.to.plot.index,], footprints.df = footprints.df,
  #                             main = paste(names(read.origin[i]), " | ", method, "| k=", k, "| ", length(reads.to.plot), " reads"))
  #       }
  #       }) -> pam.heatmaps
  #     names(pam.heatmaps) = names(read.origin)
  #   }
  # } else {
  pam.heatmaps = NULL
  #}
  
  if (plot.binary.heatmap){
    sil.res[,] %>% as_tibble() %>% mutate(indexes = 1:nrow(.), cluster = factor(cluster, levels = arrange.clusters.by.footprints(footprints.df = footprints.df))) %>% 
      arrange(cluster, desc(sil_width)) -> sorted.silhouette.results
    if (is.null(read.origin)){
      # SingleMoleculeFootprinting::PlotSM(MethSM = list(mat = binary.matrix), RegionOfInterest = RegionOfInterest, SortedReads = list(names(pam.res[rev(sorted.silhouette.results$indexes)])), sorting.strategy = "custom") -> SM.plot
      # results$SM.plot = SM.plot
    } else {
      lapply(seq_along(read.origin), function(i){
        read.indexes = which(rownames(binary.matrix) %in% read.origin[[i]])
        reads = names(pam.res[rev(sorted.silhouette.results$indexes)])[names(pam.res[rev(sorted.silhouette.results$indexes)]) %in% read.origin[[i]]]
        if(length(reads)>0){
          MethSM = list()
          MethSM[[names(read.origin)[i]]] = binary.matrix[read.indexes,]
          SortedReads = list()
          SortedReads[[names(read.origin)[i]]] = reads
          SingleMoleculeFootprinting::PlotSM(MethSM = MethSM, RegionOfInterest = RegionOfInterest_ext, 
                                             SortedReads = SortedReads, sorting.strategy = "custom")
        }
      }) -> SM.plot
      names(SM.plot) = names(read.origin)
    }
  } else {SM.plot = NULL}
  
  if(return.data){
    return(list(SW.mat = NULL, pam.res = pam.res, sil.res = NULL, footprints.df = footprints.df, TFBSs.overlapping.locus = TFBSs.overlapping.locus, footprint.plot = footprint.plot, pam.heatmaps = pam.heatmaps, SM.plot = SM.plot))
  } else {
    return(list(SW.mat = NULL, pam.res = NULL, sil.res = NULL, footprints.df = NULL, TFBSs.overlapping.locus = NULL, footprint.plot = footprint.plot, pam.heatmaps = pam.heatmaps, SM.plot = SM.plot))
  }
  
}

#' One wrapper function to go from sampleSheet + RegionOfInterest to allelic imbalance according to pam clustering
#' 
#' @param sampleSheet
#' @param RegionOfInterest will be resized in place
#' @param TFBSs will be subset for TFBSs overlapping RegionOfInterest (before resize)
#' @param params list of parameters to pass to lower functions
#' @param deduplicate TRUE or FALSE (default). Useful with samples with high PCR duplication
#' @param CytosinesToMask
#' @param SMF.plot logical
#' @param SMF.plot.params
#' 
#' @import SingleMoleculeFootprinting
#' @import magrittr
#' 
compute.pam.clustering.Wrapper = function(
    sampleSheet, Samples, genome, RegionOfInterest, TFBSs, params,
    deduplicate = FALSE, CytosinesToMask = NULL, pool.samples = TRUE, SMF.plot = FALSE){
  
  # Unpack parameters ####
  CallContextMethylation.coverage = params$CallContextMethylation[["coverage"]]
  CallContextMethylation.ConvRate.thr = params$CallContextMethylation[["ConvRate.thr"]]
  # SortReadsByTFCluster.coverage = params$SortReadsByTFCluster[["coverage"]]
  PlottingSNPs = params$PlotSingleSiteSMF[["PlottingSNPs"]]
  output.params = params$output.params
  ########################
  
  RegionOfInterest_ext = resize(RegionOfInterest, width = 500, fix = "center")
  if(!is.null(TFBSs)){
    TFBSs.overlapping.locus = plyranges::filter_by_overlaps(TFBSs[seqnames(TFBSs) == seqnames(RegionOfInterest)], RegionOfInterest_ext)
  } else {
    TFBSs.overlapping.locus = NULL
  }
  Experiment = SingleMoleculeFootprinting::DetectExperimentType(Samples = Samples)
  
  CallContextMethylation(sampleSheet = sampleSheet, sample = Samples,
                         genome = genome, RegionOfInterest = RegionOfInterest_ext,
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
  
  if(deduplicate){
    message("removing PCR duplicates and incomplete molecules")
    # source("/g/krebs/barzaghi/Rscripts/R_package/SingleMoleculeFootprinting/R/utilities.r")
    lapply(Methylation[[2]], function(x){
      incomplete.molecules.idx = as.logical(rowSums(x == 0) > 0)
      x = x[!incomplete.molecules.idx,]
      duplicated.molecules.idx = as.logical(duplicated(as.matrix(x)))
      x = x[!duplicated.molecules.idx,]
      return(x)
    }) -> Methylation[[2]]
    .MethGR = Methylation[[1]]
    Methylation[[1]] = SingleMoleculeFootprinting::MethSM.to.MethGR(MethSM =  Methylation[[2]], chromosome = as.character(unique(seqnames(RegionOfInterest_ext))))
    overlaps = GenomicRanges::findOverlaps(.MethGR, Methylation[[1]], ignore.strand = TRUE)
    strand(Methylation[[1]][subjectHits(overlaps)]) = strand(.MethGR[queryHits(overlaps)])
    Methylation[[1]]$GenomicContext = NA
    Methylation[[1]][subjectHits(overlaps)]$GenomicContext = .MethGR[queryHits(overlaps)]$GenomicContext
    colnames(elementMetadata(Methylation[[1]])) = gsub("_T$", "_Coverage", colnames(elementMetadata(Methylation[[1]])))
    colnames(elementMetadata(Methylation[[1]])) = gsub("_M$", "_MethRate", colnames(elementMetadata(Methylation[[1]])))
    elementMetadata(Methylation[[1]])[grep("_MethRate$", colnames(elementMetadata(Methylation[[1]])))] =
      as.matrix(elementMetadata(Methylation[[1]])[grep("_MethRate$", colnames(elementMetadata(Methylation[[1]])))]) /
      as.matrix(elementMetadata(Methylation[[1]])[grep("_Coverage$", colnames(elementMetadata(Methylation[[1]])))])
  }
  
  if (length(Methylation[[1]]) == 0){
    message("No cytosines found covering the region")
    return()
  }
  
  if (!is.null(PlottingSNPs)){
    PlottingSNPs = PlottingSNPs[queryHits(findOverlaps(PlottingSNPs, GRanges(unique(seqnames(Methylation[[1]])), IRanges(min(start(Methylation[[1]])), max(end(Methylation[[1]]))))))]
    params$PlotSingleSiteSMF[["PlottingSNPs"]] = PlottingSNPs
  }
  
  if (SMF.plot){
    PlotAvgSMF(MethGR = Methylation[[1]], RegionOfInterest = RegionOfInterest_ext, ShowContext = TRUE,
               TFBSs = TFBSs.overlapping.locus, SNPs = PlottingSNPs) -> SingleSite.plot
  } else {SingleSite.plot = NULL}
  
  read.origin = lapply(Methylation[[2]], rownames)
  binary.matrix = Reduce(SingleMoleculeFootprinting::rbind.fill.matrix.sparse, Methylation[[2]])
  compute.pam.clustering(binary.matrix = binary.matrix, TFBSs.overlapping.locus = TFBSs.overlapping.locus, RegionOfInterest = RegionOfInterest, RegionOfInterest_ext = RegionOfInterest_ext,
                         params = params, read.origin = read.origin,
                         output.params = output.params) -> final.results

  Reduce(rbind,
  lapply(seq_along(read.origin), function(i){

    current.sample = names(read.origin)[i]
    current.reads = read.origin[[i]]
    all.cluster.nrs = seq(params$pam[["k"]])
    clustering.reads = final.results$pam.res[names(final.results$pam.res) %in% current.reads]
    
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

  })) -> cluster.coverage.df

  final.results$cluster.coverage.df = cluster.coverage.df
  final.results$read.origin = read.origin
  final.results$SingleSite.plot = SingleSite.plot
  final.results$RegionOfInterest = RegionOfInterest
  
  return(final.results)
  
}
