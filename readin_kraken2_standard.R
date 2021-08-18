readin_kraken2_standard <- function(directory, recursive = T, saveIntermediate = F){
  
  require(data.table)
  require(tidyverse)
  
  print(paste0("Looking for kraken2 standard output to load..."))
  
  files <- Sys.glob(file.path(directory, "*standard.txt"))
  
  dflist <- NULL
  
  print(paste0("Found ", length(files), " files to load..."))
  
  for(i in 1:length(files)) {
    print(paste0("  >> loading ... ", files[i]))
    
    dflist[[i]] <- fread(files[i], header = F,
                         col.names = c("uc","read","taxon","length","LCAmap")) %>%
      mutate(filename=files[i]) %>%
      dplyr::count(filename, uc, taxon, name = "numseqs", sort = T)
    
    if (saveIntermediate){
      saveRDS(dflist[[i]], paste0(files[i],".rds"))
     
    }
  }
  
  kraken2 <- do.call(rbind,dflist) %>%
    mutate(filename = sapply(str_split(filename,"/"), function(x) x[length(x)]),
           seq_id=gsub("\\.standard\\.txt","",filename),
           taxidex = str_extract(taxon,"taxid [0-9]+")) %>%
    separate(taxidex, c("taxfoo","taxid"), convert = T) %>%
    select(seq_id, uc, taxon, taxid, numseqs, filename)
  
  return(kraken2)
}