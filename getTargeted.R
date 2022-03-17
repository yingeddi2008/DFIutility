argopts <- commandArgs(trailingOnly = TRUE)
path = argopts[1]

if (is.na(argopts[2])){
   sidx = 8
  } else {
   sidx = as.numeric(argopts[2])
}

if (is.na(argopts[3])){
   intermediate = F
   cat("No intermediate files will be outputted!\n")
  } else {
   intermediate = as.logical(argopts[3])
   if (intermediate) cat("Intermediate files will be outputted!\n")
}

# cat(intermediate)

library(data.table)
library(tidyverse)

try(source("~/Documents/Eddi/DFIutility/readDiamond.R"))
try(source("~/OneDrive - The University of Chicago/DFIutility/readDiamond.R"))


dmndfns <- dir(path = path, pattern="dmnd", recursive = T, full.names = T)
nfn <- length(dmndfns)
cat("Found ", nfn," files to parse!\n")

# dmndfns <- Sys.glob(path)

tlist = NULL

for (i in 1:nfn){
  d = dmndfns[i]
  
  if (file.size(d) < 300) { next } 
  
  cat(paste0("Progress: ",i,"/",nfn, ". "))
  
  tmp <- readDiamond(d)
  
  if (is.null(tmp)) { next }

  ftmp <- tmp %>%
    ungroup() %>%
    dplyr::count(filename, sseqid, slen, name = "rawCnt") %>%
    mutate(seq_id = basename(filename),
           seq_id = gsub("\\.R[12].[a-zA-Z_]+.dmnd","", seq_id),
           dbsource = sapply(str_split(filename,"/"), function(x) x[sidx])) %>%
    select(seq_id, dbsource, sseqid, slen, rawCnt, filename)
  
  tlist[[d]] <- ftmp

if (intermediate){  
  saveRDS(ftmp, paste0(basename(d),".rds"))
  }
}

big <- dplyr::bind_rows(tlist) %>%
  dplyr::count(seq_id, dbsource, sseqid, slen, wt = rawCnt, name = "mappedReads")

saveRDS(big,"targeted.dnmd.rds")
