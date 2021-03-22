argopts <- commandArgs(trailingOnly = TRUE) 
path = argopts[1]

try(source("~/Documents/Eddi/DFIutility/readin_emapper.R"))
try(source("~/OneDrive - The University of Chicago/DFIutility/readin_emapper.R"))
kall = readin_emapper(path, recursive = F)
saveRDS(kall, "all.emapper.rds")
