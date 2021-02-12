readin_emapper <- function(directory, recursive = F,
                           pattern="emapper.annotations$" ){
  # a read-in function for emapper annotations output
  require(tidyverse)
  
  cat("Looking for emapper annotation files to load...\n")
  
  efns <- dir(path = directory, recursive = recursive, 
              pattern=pattern,
              full.names = TRUE)
  
  cat("Found", length(efns), "files to load...\n")
  
  eflist <- list()
  
  for (ef in efns){
    cat("  >> Loading ", ef, " ...\n")
    eflist[[ef]] <- read_tsv(ef, comment = "#", 
                             col_names = c("query_name","seed_eggNOG_ortholog", "seed_ortholog_evalue",
                                           "seed_ortholog_score","best_tax_level","Preferred_name",
                                           "GOs","EC","KEGG_ko","KEGG_Pathway","KEGG_Module","KEGG_Reaction",
                                           "KEGG_rclass","BRITE","KEGG_TC",
                                           "CAZy","BiGG_Reaction","taxonomic scope",
                                           "eggNOG OGs","best eggNOG OG","COG Functional cat.",
                                           "eggNOG free text desc."),
                             col_types = "ccnncccccccccccccccccc") %>%
      mutate(sample = gsub(pattern = pattern,"",basename(ef)))
  }
  
  efall <- bind_rows(eflist)
  colnames(efall)[1] <- "query_name"
  
  return(efall)
}