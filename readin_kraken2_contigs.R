#readin kraken2 function
readin_kraken2_contigs <- function(directory){
  
  files <- dir(directory, pattern="standard")
  dflist <- NULL
  for(i in 1:length(files)) {
    print(i)
    int <- fread(file.path(directory,files[i])) %>%
      mutate(filename=files[i])
    
    dflist[[i]] <- int
  }
  kraken2 <- do.call(rbind,dflist) %>%
    select(contig=V2,taxon=V3,length=V4,filename) %>%
    mutate(seq_id=gsub("\\.standard\\.txt","",filename),
           cohort="biobank")
  return(kraken2)
}
