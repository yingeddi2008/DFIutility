readin_json_crisprFinder <- function(directory="."){
  
  js <- Sys.glob(paste0(directory,"/*/result.json"))
  
  require(jsonlite)
  
  cat(paste0("There are ",length(js), " JSON files to parse ...\n"))
  
  pList <- list()
  
  for (f in js){
    
    cat(paste0("Processing ",f," ...\n"))
    isoName <- str_extract(f,"ST[0-9]+_[0-9]+")
    
    if (file.size(f) < 4500) next()
    
    parsedCasInt <- fromJSON(f)[[4]] %>%
      as_tibble() %>%
      unnest(Cas, keep_empty = T) 
    
    # check if cas-system was found
    if (all(c("Start","End","Type","Genes") %in% names(parsedCasInt))) {
      
      parsedInt <- parsedCasInt %>% # view()
        dplyr::rename(casStart = Start,
                      casEnd = End,
                      casType = Type)  %>% # view
        select(-Genes) %>%  #### do not expand genes
        unnest(Crisprs) # removed any contigs that do not have crisprs
      
    } else {
      
      parsedInt <- parsedCasInt %>% # view()
        dplyr::mutate(casStart = NA,
                      casEnd = NA,
                      casType = NA)  %>% # view
        unnest(Crisprs)
      
    }
    
    
    if (nrow(parsedInt) == 0) next() #
    
    parsed <- parsedInt %>%   
      dplyr::rename(crisprStart = Start, 
                    crisprEnd = End,
                    crisprAT = AT,
      ) %>% # view
      unnest(Regions) %>%   
      dplyr::rename(regionStart = Start, 
                    regionEnd = End,
                    regionAT = AT,
                    regionType = Type,
                    regionSequence = Sequence) %>% #view
      mutate(isolate = isoName)
    
    # expand genes if want to
    # mutate(Genes = paste0(Genes, collapse = ";")) %>%
    # select(Genes)
    # unnest(Genes, keep_empty = T) %>% # view
    # dplyr::rename(geneStart = Start,
    #               geneEnd = End) %>%
    # mutate(isolate = isoName) %>%
    # select(isolate, contigName = Name, starts_with("cas"),
    #        starts_with("crispr"), starts_with("region"), everything(.))
    
    pList[[isoName]] <- parsed
  }
  
  return(bind_rows(pList))
  
}
