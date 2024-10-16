library(SingleMoleculeFootprinting)
library(BSgenome.Mmusculus.UCSC.mm10)
library(testthat)

CacheDir <- ExperimentHub::getExperimentHubOption(arg = "CACHE")
Qinput = paste0(CacheDir, "/NRF1Pair_sampleFile.txt")

test_that("BaitCapture returns a numeric", {

  BaitCapture(
    sampleFile = Qinput,
    genome = BSgenome.Mmusculus.UCSC.mm10,
    baits = SingleMoleculeFootprintingData::EnrichmentRegions_mm10.rds()) %>%
    is.numeric() %>%
    expect_true()

})
