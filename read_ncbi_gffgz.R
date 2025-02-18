readin_ncbi_gffgz <- function(directory, recursive=T){
  require(tidyverse)
  require(reshape2)
  require(rtracklayer)
  require(stringr)
  require(yingtools2)
  require(seqinr)
  require(dplyr)
  
  print(paste0("Looking for GFF files to load..."))
  
  gffs <- Sys.glob(directory)
  dflist <- NULL
  col_list <- NULL

  print(paste0("Found ", length(gffs), " files to load..."))
  
  for(i in 1:length(gffs)){
    print(paste0("loading ... ", gffs[i]))
    #i=1
    
    df <- readGFF(gzfile(gffs[i])) %>% 
      as.data.frame() %>%
      dplyr::select(any_of(c("seqid","source","type","start", 
                             "end", "score","strand", "phase",
                             "locus_tag","product","note",
                             "db_xref","gene", "eC_number"))) %>%
      mutate(filename=gffs[i]) 
    
    df2 <- cbind(i,ncol(df))
    
    dflist[[i]] <- df
    col_list[[i]] <- df2
  }
  
  big_gff <- bind_rows(dflist) 
  
  return(big_gff)
  
}
