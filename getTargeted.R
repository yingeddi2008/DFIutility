argopts <- commandArgs(trailingOnly = TRUE)
path = argopts[1]
sidx = as.numeric(argopts[2])

if (is.null(argopts[3])){
   intermediate = F
  } 
else {
    intermediate = as.logic(argopts[3])
}

library(data.table)
library(tidyverse)

source("~/Documents/Eddi/funcsNscripts/readDiamond.R")

dmndfns <- Sys.glob(path)

tlist = NULL

for (d in dmndfns){
  
  if (file.size(d) < 300) { next }   
  tmp <- readDiamond(d)
  
  if (is.null(tmp)) { next }

  ftmp <- tmp %>%
    ungroup() %>%
    dplyr::count(filename, sseqid, slen, name = "rawCnt") %>%
    mutate(seq_id = basename(filename),
           seq_id = sapply(str_split(seq_id,"_"), function(x) x[1]),
           dbsource = sapply(str_split(filename,"/"), function(x) x[sidx])) %>%
    select(seq_id, dbsource, sseqid, slen, rawCnt, filename)
  
  tlist[[d]] <- ftmp

if (intermediate){  
  saveRDS(ftmp, paste0(basename(d),".rds"))
  }
}

big <- bind_rows(tlist) %>%
  count(seq_id, dbsource, sseqid, slen, wt = rawCnt, name = "mappedReads")

saveRDS(big,"targeted.dnmd.rds")
