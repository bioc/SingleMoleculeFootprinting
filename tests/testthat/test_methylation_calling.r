library(SingleMoleculeFootprinting)
library(BSgenome.Mmusculus.UCSC.mm10)
library(testthat)

CacheDir <- ExperimentHub::getExperimentHubOption(arg = "CACHE")
Qinput = paste0(CacheDir, "/NRF1Pair_sampleFile.txt")
MySample <- suppressMessages(readr::read_delim(Qinput, delim = "\t")[[2]])
RegionOfInterest <- GRanges("chr6", IRanges(start = 88106000, end = 88106500))
QuasRprj = GetQuasRprj(Qinput, BSgenome.Mmusculus.UCSC.mm10)

test_that("GetSingleMolMethMat returns a matrix", {

  GetSingleMolMethMat(QuasRprj, range = RegionOfInterest, sample = MySample)$NRF1pair_DE_ %>%
    inherits(., "dgCMatrix") %>%
    expect_true()

})

test_that("DetectExperimentType behaves as expected", {

  MySample = QuasR::alignments(QuasRprj)[[1]]$SampleName
  DetectExperimentType(Samples = MySample) %>%
    is.character() %>%
    expect_true()
  DetectExperimentType(Samples = MySample) %>%
    expect_message()

})


