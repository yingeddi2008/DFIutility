read_tabBlastn <- function(filename, ident = 85, cov = 85 ){
  
  require(tidyverse)
  require(readr)
  require(dplyr)
  
  CNAMES <- c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend", 
             "sstart", "send","evalue", "bitscore", "score", "qlen", "qcovus", "sstrand", "slen")
  
  bo <- readr::read_tsv(filename, col_names = F) 
  
  colnames(bo) <- CNAMES[1:ncol(bo)]
  
  bo %>%
    mutate(qcov = (qend - qstart + 1)/qlen * 100,
           iterPos = ifelse(sstrand == "minus", send, sstart),   #switch start and end position in minus
           send = ifelse(sstrand == "minus", sstart, send),
           sstart = ifelse(sstrand == "minus", iterPos, sstart),
           iterPos = NULL,
           filename = filename) %>%
    dplyr::filter(pident > ident, qcov > cov)
  
}