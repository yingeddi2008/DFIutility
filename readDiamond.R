readDiamond <- function(diaFn, ident = 80, ecut = 1e-5, qlenP = 0.8){
  require(data.table)
  require(tidyverse)

  cat("loading ..", diaFn, "\n")
  
  if (file.size(diaFn) < 500 ){
    warning("Skipping ", diaFn, " due to small file size!")
    return()
  }
  
  tmp <- fread(diaFn, header = F, 
               col.names = c("qseqid","sseqid","pident","length","mismatch",
                             "gapoen","qstart","qend","sstart","send",
                             "evalue","bitscore","qlen","slen")) %>%
    mutate(filename=diaFn)

  # 80pct identity, 80% translated read length, best hit per read
  filint <- tmp %>%
    filter(pident >= ident,
           evalue < ecut,
           length > qlenP*qlen/3) 
  
  if (nrow(filint) > 0) {
    
    filint %>%
      group_by(qseqid) %>%
      top_n(1, wt = bitscore) %>%
      top_n(1, wt = pident) %>%
      top_n(-1, wt = evalue) %>%
      add_count(sseqid, name = "swt") %>%
      sample_n(1, weight = swt)
  }
  
  
}
