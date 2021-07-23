readin_prokka_tables <- function(directory){
  #readin prokka annotation and seq tables
  #returns list with two dataframes, 1 is annotation, 2 is sequence with locus tag
  require(tidyverse)
  require(reshape2)
  require(rtracklayer)
  require(stringr)
  require(yingtools2)
  require(seqinr)
  
  print(paste0("Looking for prokka annotation files to load..."))
  
  gene_files <- Sys.glob(file.path(directory,"*/*tsv"))
  nuc_seq <- Sys.glob(file.path(directory,"*/*ffn"))
  prot_seq <- Sys.glob(file.path(directory,"*/*faa"))
  
  print(paste0("Found ", length(nuc_seq), " files to load..."))
  
  #make gene table & seq table
  dflist <- NULL
  seqlist <- NULL
  for(i in 1:length(gene_files)){
    print(paste0("loading ... ", file.path(gene_files[i])))
    #i=1
    int <- read.delim(file=file.path( gene_files[i]),stringsAsFactors = F) %>%
      mutate(filename=gene_files[i])
    #glimpse(int)
    
    dflist[[i]] <- int
    
    print(paste0("loading ... ", file.path(nuc_seq[i])))
    intn <- read.delim(file=file.path(nuc_seq[i]),stringsAsFactors = F,header=F) %>%
      mutate(seq_name=str_extract(V1,">.+")) %>%
      fill(seq_name) %>%
      filter(!grepl(V1,pattern="^>")) %>%
      group_by(seq_name) %>%
      summarize(nuc_sequence=paste0(V1,collapse="")) %>%
      ungroup() %>%
      mutate(locus_tag=str_extract(seq_name,pattern="^>[A-Z]+\\_[0-9]+"),
             locus_tag=gsub(">","",locus_tag),
             product=gsub("^>[A-Z]+\\_[0-9]+ ","",seq_name),
             filename=gsub(".ffn","",nuc_seq[i]))
    
    intp <- read.delim(file=file.path(prot_seq[i]), stringsAsFactors = F,header=F) %>%
      mutate(seq_name=str_extract(V1,">.+")) %>%
      fill(seq_name) %>%
      filter(!grepl(V1,pattern="^>")) %>%
      group_by(seq_name) %>%
      summarize(prot_sequence=paste0(V1,collapse="")) %>%
      ungroup() %>%
      mutate(locus_tag=str_extract(seq_name,pattern="^>[A-Z]+\\_[0-9]+"),
             locus_tag=gsub(">","",locus_tag),
             product=gsub("^>[A-Z]+\\_[0-9]+ ","",seq_name),
             filename=gsub(".faa","",prot_seq[i]))
    
    intseq <- intn %>%
      left_join(intp)
    
    seqlist[[i]] <- intseq
    
  }
  
  prokka <- do.call(rbind,dflist)
  prokka_seq <- do.call(rbind,seqlist)
  
  dflist <- NULL
  dflist[[1]] <- prokka
  dflist[[2]] <- prokka_seq
  return(dflist)
}