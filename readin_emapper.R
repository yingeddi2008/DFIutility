readin_emapper <- function(directory, recursive = F){
  # a read-in function for emapper annotations output
  require(tidyverse)
  
  cat("Looking for emapper annotation files to load...\n")
  
  efns <- dir(path = directory, recursive = recursive, 
              pattern="emapper.annotations$",
              full.names = TRUE)
  
  cat("Found", length(efns), "files to load...\n")
  
  eflist <- list()
  
  for (ef in efns){
    eflist[[ef]] <- read_tsv(ef, skip = 3) %>%
      mutate(sample = gsub("emapper.annotations","",basename(ef)))
  }
  
  efall <- bind_rows(eflist)
  colnames(efall)[1] <- "query_name"
  
  return(efall)
}