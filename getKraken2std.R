argopts <- commandArgs(trailingOnly = TRUE) 
path = argopts[1]
saveF = as.logical(argopts[2])

try(source("~/Documents/Eddi/DFIutility/readin_kraken2_standard.R"))
try(source("~/OneDrive - The University of Chicago/DFIutility/readin_kraken2_standard.R"))
kall = readin_kraken2_standard(path, saveIntermediate = saveF)
saveRDS(kall, "all.kraken2.standard.rds")
