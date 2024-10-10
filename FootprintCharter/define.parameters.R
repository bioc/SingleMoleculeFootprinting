#' #### Convenience stuff that I don't know where else to put
#' #' 
#' #' @param tax_group defaults to "vertebrates"
#' #' 
#' #' @import JASPAR2020
#' #' @import TFBSTools
#' library(JASPAR2020)
#' library(TFBSTools)
#' prepare.JASPAR.PWMs = function(tax_group="vertebrates"){
#'   
#'   opts <- list()
#'   opts[["tax_group"]] <- tax_group
#'   PFMatrixList <- getMatrixSet(JASPAR2020, opts)
#'   PWMatrixList <- toPWM(PFMatrixList, pseudocounts=0.8)
#'   
#'   return(PWMatrixList)
#'   
#' }
#' 
PWM.list = NULL#prepare.JASPAR.PWMs(tax_group = "vertebrates")

archetypes = NULL#Load.JASPAR.archetypes()

##################################
# Amplicon parameters
##################################

amplicon.parameters.list = list(
  filter.single.molecule.matrix = c(NA.thr.rows = 0, NA.thr.cols = 0),
  make.full.coordinate.matrix = c(padding = 15), # should be half of window.size
  matrix.sliding.window.average = c(window.size=40),
  fill.empty.columns = c(max.nr.cols = 30),
  unsupervised.clustering.coverage = 20, # at the very least, should be higher than pam[["k"]]
  compute.pam.clustering = c(max.nr.reads = list(NULL), detect.fp.from.binary = TRUE), # set to NULL will allow all the reads to be used
  parDist = c(method = "euclidean", threads = 1),
  pam = c(k = 24),
  detect.footprint = c(TF.fp.min.length = 5, TF.fp.max.length = 75, Nucleosome.fp.min.length = 120, Nucleosome.fp.max.length = 1000, cytosine.min.coverage = 5), # TF.fp.min.length is the shortest JASPAR2020 PWM
  # arrange.clusters.by.footprints = c(min.fp.thr = 0.25, TF.fp.min.length = 15, TF.fp.max.length = 50, Nucleosome.fp.min.length = 140), # 147
  annotate.orphan.footprints = c(PWM.list=PWM.list, archetypes=list(archetypes), genome=NULL, use.passed.TFBSs.for.annotation=TRUE, scan.orphan.footprints=FALSE, cores = 1), # cores here is an unsolved problem: when using the multisite wrapper this will ~multiply by the number of cores passed there
  CallContextMethylation = c(coverage = 40, ConvRate.thr = list(NULL)),
  SortReadsByTFCluster = c(coverage = 30),
  PlotSingleSiteSMF = c(PlottingSNPs = list(NULL)), # GRanges object of SNPs to visualize. Assumed to be already subset. Assumed to have the reference and alternative sequences respectively under the columns "R" and "A"
  misc = c(nanopore.pairs.filtering = FALSE),
  output.params = list(plot.SW.heatmap = FALSE, plot.binary.heatmap = FALSE, compute.silhouette = FALSE, return.data = TRUE)
)

##################################
# Bait capture specific parameters
##################################

bait.capture.parameters.list = list(
  filter.single.molecule.matrix = c(NA.thr.rows = 0, NA.thr.cols = 0),
  make.full.coordinate.matrix = c(padding = 20), # should be half of window.size
  matrix.sliding.window.average = c(window.size=40),
  fill.empty.columns = c(max.nr.cols = 20),
  unsupervised.clustering.coverage = 20, # at the very least, should be higher than pam[{"k"}]
  compute.pam.clustering = c(max.nr.reads = list(NULL), detect.fp.from.binary = TRUE), # set to NULL will allow all the reads to be used
  parDist = c(method = "euclidean", threads = 1),
  pam = c(k = 12),
  detect.footprint = c(TF.fp.min.length = 5, TF.fp.max.length = 75, Nucleosome.fp.min.length = 120, Nucleosome.fp.max.length = 1000, cytosine.min.coverage = 5), # TF.fp.min.length is the shortest JASPAR2020 PWM
  # arrange.clusters.by.footprints = c(min.fp.thr = 0.25, TF.fp.min.length = 15, TF.fp.max.length = 50, Nucleosome.fp.min.length = 140), # 147
  annotate.orphan.footprints = c(PWM.list=PWM.list, archetypes=list(archetypes), genome=NULL, use.passed.TFBSs.for.annotation=TRUE, scan.orphan.footprints=FALSE, cores = 1),
  CallContextMethylation = c(coverage = 20, ConvRate.thr = list(NULL)),
  SortReadsByTFCluster = c(coverage=20),
  PlotSingleSiteSMF = c(PlottingSNPs = list(NULL)), # GRanges object of SNPs to visualize. Assumed to be already subset. Assumed to have the reference and alternative sequences respectively under the columns "R" and "A"
  misc = c(nanopore.pairs.filtering = FALSE),
  output.params = list(plot.SW.heatmap = FALSE, plot.binary.heatmap = FALSE, compute.silhouette = FALSE, return.data = TRUE)
)

##################################
# Nanopore specific parameters: haven't checked them in a long time, be very careful if using
##################################

# nanopore.parameters.list = list(
#   filter.single.molecule.matrix = c(NA.thr.rows = 0, NA.thr.cols = 0),
#   make.full.coordinate.matrix = c(padding = 20), # should be half of window.size
#   matrix.sliding.window.average = c(window.size=40),
#   fill.empty.columns = c(max.nr.cols = 20),
#   unsupervised.clustering.coverage = 20, # at the very least, should be higher than pam[{"k"}]
#   compute.pam.clustering = c(max.nr.reads = list(NULL), detect.fp.from.binary = TRUE), # set to NULL will allow all the reads to be used
#   parDist = c(method = "euclidean", threads = 1),
#   pam = c(k = 10),
#   detect.footprint = c(TF.fp.min.length = 5, TF.fp.max.length = 75, Nucleosome.fp.min.length = 120, cytosine.min.coverage = 5), # TF.fp.min.length is the shortest JASPAR2020 PWM
#   # arrange.clusters.by.footprints = c(min.fp.thr = 0.25, TF.fp.min.length = 15, TF.fp.max.length = 50, Nucleosome.fp.min.length = 140), # 147
#   annotate.orphan.footprints = c(PWM.list=PWM.list, archetypes=list(archetypes), genome=BSgenome.Dmelanogaster.UCSC.dm6, use.passed.TFBSs.for.annotation=TRUE, scan.orphan.footprints=FALSE, cores = 1),
#   CallContextMethylation = c(coverage = 20, ConvRate.thr = list(NULL)),
#   SortReadsByTFCluster = c(coverage=20),
#   PlotSingleSiteSMF = c(PlottingSNPs = list(NULL)), # GRanges object of SNPs to visualize. Assumed to be already subset. Assumed to have the reference and alternative sequences respectively under the columns "R" and "A"
#   misc = c(nanopore.pairs.filtering = FALSE),
#   output.params = list(plot.SW.heatmap = FALSE, plot.binary.heatmap = FALSE, compute.silhouette = TRUE, return.data = TRUE)
# )



