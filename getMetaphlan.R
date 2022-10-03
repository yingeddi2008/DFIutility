argopts <- commandArgs(trailingOnly = TRUE)
path = argopts[1]

library(data.table)
library(tidyverse)

readin_bracken_reports <- function(directory){
  files <- dir(directory, pattern="metagenome.txt")
  
  dflist <- NULL
  
  for(i in 1:length(files)) {
    int <- read_tsv(file.path(directory,files[i]), skip = 5) %>%
      mutate(filename = files[i]) %>%
      mutate(seq_id = gsub("_metagenome.txt", "", filename),
             taxid = sapply(str_split(clade_taxid, "\\|"), function(x) x[length(x)]))

    dflist[[i]] <- int
  }
  print(paste("Total files read:",i))
  mpa <- bind_rows(dflist) 
  return(mpa)
}

kall = readin_bracken_reports(path)
saveRDS(kall, "all.metaphlan.rds")
