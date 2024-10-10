library(SingleMoleculeFootprinting)
library(BSgenome.Mmusculus.UCSC.mm10, lib.loc = "/g/krebs/barzaghi/R/x86_64-pc-linux-gnu-library/4.1/")
library(tidyverse)
library(magrittr)
library(parallel)

source("/g/krebs/barzaghi/Rscripts/CrappyUtils.R")
source("/g/krebs/barzaghi/analyses/single.molecule.classification/utils/define.parameters.R") # <- For Mathias
source("/g/krebs/barzaghi/analyses/single.molecule.classification/utils/create.sliding.windows.R") # <- For Mathias
source("/g/krebs/barzaghi/analyses/single.molecule.classification/utils/detect.footprints.R")
source("/g/krebs/barzaghi/analyses/single.molecule.classification/utils/annotate.footprints.R")
source("/g/krebs/barzaghi/analyses/single.molecule.classification/utils/pam.clustering.R") # <- For Mathias (only compute.pam.clustering function)
source("/g/krebs/barzaghi/analyses/single.molecule.classification/utils/MultiSite.wrapper.R")
source("/g/krebs/barzaghi/analyses/single.molecule.classification/utils/master.table.R")
source("/g/krebs/barzaghi/analyses/single.molecule.classification/utils/master.table_utils.R")
source("/g/krebs/barzaghi/analyses/single.molecule.classification/utils/patchwork.single.site.R")
