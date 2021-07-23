readin_prokka_gff <- function(directory){
  require(tidyverse)
  require(reshape2)
  require(rtracklayer)
  require(stringr)
  require(yingtools2)
  require(seqinr)
  require(dplyr)
  
  print(paste0("Looking for GFF files to load..."))
  
  gffs <- Sys.glob(file.path(directory,"*/*gff"))
  dflist <- NULL
  col_list <- NULL
  
  print(paste0("Found ", length(gffs), " files to load..."))
  
  for(i in 1:length(gffs)){
    print(paste0("loading ... ", file.path(gffs[i])))
    #i=1
    df <- readGFF(file.path(gffs[i])) %>% 
      as.data.frame() %>%
      dplyr::select(any_of(c("seqid","source","type","start", 
                      "end", "score","strand", "phase",
                      "locus_tag","product","note",
                      "db_xref","gene", "eC_number","ID","Name"))) %>%
      mutate(filename=gffs[i])
    
    df2 <- cbind(i,ncol(df))
    
    dflist[[i]] <- df
    col_list[[i]] <- df2
  }
  
  big_gff <- bind_rows(dflist) 
  
  # inf <- do.call(rbind,dflist) %>%
  #   mutate(row=row_number()) %>%
  #   group_by(row) %>%
  #   summarize(inf=paste(unlist(inference),collapse="|"))
  # 
  # big_gff_final <- big_gff %>%
  #   left_join(inf) %>%
  #   mutate(uni_full=str_extract(inf,pattern="UniProtKB.+"),
  #          uni_num=gsub("UniProtKB\\:","",uni_full))
  
  return(big_gff)
  # return(big_gff_final)
  
}
