argopts <- commandArgs(trailingOnly = TRUE) 
path = argopts[1]

source("~/Documents/Eddi/DFIutility/readin_kraken2_standard.R")
kall = readin_kraken2_standard(path, recursive = F, saveIntermediate = T)
saveRDS(kall, "all.kraken2.standard.rds")
