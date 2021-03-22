argopts <- commandArgs(trailingOnly = TRUE) 
path = argopts[1]

try(source("~/Documents/Eddi/DFIutility/readin_kraken2_standard.R"))
try(source("~/OneDrive - The University of Chicago/DFIutility/readin_kraken2_standard.R"))
kall = readin_kraken2_standard(path, recursive = F, saveIntermediate = F)
saveRDS(kall, "all.kraken2.standard.rds")
