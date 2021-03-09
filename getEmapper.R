argopts <- commandArgs(trailingOnly = TRUE) 
path = argopts[1]

source("~/Documents/Eddi/DFIutility/readin_emapper.R")
kall = readin_emapper(path, recursive = F)
saveRDS(kall, "all.emapper.rds")
