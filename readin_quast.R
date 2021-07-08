readin_quast <- function(directory, recursive=T){
  # readin 
  
  require(tidyverse)
  
  print(paste0("Looking for quast tsv files to load..."))
  
  qtsvs <- dir(path = directory, recursive = recursive, pattern="transposed_report.tsv")

  print(paste0("Found ", length(qtsvs), " files to load..."))
  
  # place holder for quast
  aqst <- lapply(qtsvs, read_tsv)
  
  longquast <- bind_rows(aqst) %>%
    mutate(Assembly = gsub("[_\\.]spades.contigs","",Assembly),
           Assembly = gsub("_","-",Assembly)) %>%
    filter(grepl("^EP|DFI|MSK|SL",Assembly)) %>%
    gather("stat","value",-Assembly) %>%
    arrange(Assembly) %>%
    dplyr::rename(seq_id = Assembly)

  return(longquast)
}