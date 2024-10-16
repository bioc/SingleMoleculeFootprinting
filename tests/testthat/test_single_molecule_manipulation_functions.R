library(SingleMoleculeFootprinting)
library(BSgenome.Mmusculus.UCSC.mm10)
library(testthat)

CacheDir <- ExperimentHub::getExperimentHubOption(arg = "CACHE")
sampleFile = paste0(CacheDir, "/NRF1Pair_sampleFile.txt")
samples <- suppressMessages(unique(readr::read_delim(sampleFile, delim = "\t")$SampleName))
RegionOfInterest <- GRanges("chr6", IRanges(88106000, 88106500))
CallContextMethylation(
  sampleFile = sampleFile, 
  samples = samples, 
  genome = BSgenome.Mmusculus.UCSC.mm10, 
  RegionOfInterest = RegionOfInterest,
  coverage = 20, 
  returnSM = TRUE,
  ConvRate.thr = NULL,
  verbose = FALSE,
  clObj = NULL
) -> Methylation
TFBSs = qs::qread(system.file("extdata", "TFBSs_1.qs", package="SingleMoleculeFootprinting"))


test_that("BinMethylation returns named numeric vector", {

  TFBS_cluster = sort(TFBSs, by = ~ seqnames + start + end)
  TFBS_centers = start(TFBSs) + (end(TFBSs)-start(TFBSs))/2
  BinsCoordinates = IRanges(TFBS_centers-7, TFBS_centers+7)
  BinMethylation(MethSM = Methylation[[2]]$NRF1pair_DE_, Bin = BinsCoordinates) %>%
    is.numeric() %>%
    expect_true()

})

test_that("SortReads returns a list", {

  bins = list(c(-35,-25), c(-7,7), c(25,35))
  TFBS_cluster = sort(TFBSs, by = ~ seqnames + start + end)
  TFBS_centers = start(TFBS_cluster) + (end(TFBS_cluster)-start(TFBS_cluster))/2
  BinsCoordinates = IRanges(start = c(min(TFBS_centers)+bins[[1]][1], TFBS_centers+bins[[2]][1], max(TFBS_centers)+bins[[3]][1]),
                            end = c(min(TFBS_centers)+bins[[1]][2], TFBS_centers+bins[[2]][2], max(TFBS_centers)+bins[[3]][2]))
  lapply(Methylation[[2]], SortReads, BinsCoordinates = BinsCoordinates, coverage = 30) %>%
    is.list() %>%
    expect_true()

})
