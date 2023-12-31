% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/context_methylation_functions.r
\name{FilterByConversionRate}
\alias{FilterByConversionRate}
\title{Calculate reads conversion rate}
\usage{
FilterByConversionRate(MethSM, chr, genome, thr = 0.2, verbose = TRUE)
}
\arguments{
\item{MethSM}{as comes out of the func GetSingleMolMethMat}

\item{chr}{Chromosome, MethSM doesn't carry this info}

\item{genome}{BSgenome}

\item{thr}{Double between 0 and 1. Threshold above which to filter reads. Defaults to 0.2}

\item{verbose}{TRUE/FALSE}
}
\value{
Filtered MethSM
}
\description{
Calculate reads conversion rate
}
\examples{
Qinput = paste0(tempdir(), "/NRF1Pair_Qinput.txt")
library(BSgenome.Mmusculus.UCSC.mm10)

if(file.exists(Qinput)){
    QuasRprj = GetQuasRprj(Qinput, BSgenome.Mmusculus.UCSC.mm10)

    sample = readr::read_delim(Qinput, delim = "\t")$SampleName
    range = GRanges(seqnames = "chr6", ranges = IRanges(start = 88106000, end = 88106500), strand = "*")

    MethSM = GetSingleMolMethMat(QuasRprj, range, sample)
    MethSM = FilterByConversionRate(MethSM, chr = "chr6", genome = BSgenome.Mmusculus.UCSC.mm10, thr = 0.8)
}

}
