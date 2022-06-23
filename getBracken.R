argopts <- commandArgs(trailingOnly = TRUE)
path = argopts[1]

library(data.table)
library(tidyverse)

readin_bracken_reports <- function(directory){
  files <- dir(directory, pattern="report")
  dflist <- NULL
  for(i in 1:length(files)) {
    int <- fread(file.path(directory,files[i])) %>%
      mutate(filename=files[i])

    dflist[[i]] <- int
  }
  print(paste("Total files read:",i))
  bracken <- do.call(rbind,dflist) %>%
    mutate(seq_id=gsub("_report_[A-Za-z]+.txt","",filename))
  return(bracken)
}

kall = readin_bracken_reports(path, saveIntermediate = saveF)
saveRDS(kall, "all.bracken.rds")
