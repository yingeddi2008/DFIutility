getKeggPal <- function (keggtax, rn = 3, paln = 10) {
  
  require(yingtools2)
  require(tidyverse)
  
  basecols <- tibble( rank2 = c("09182 Protein families: genetic information processing",
                                "09183 Protein families: signaling and cellular processes",
                                "09181 Protein families: metabolism","09132 Signal transduction",
                                "09101 Carbohydrate metabolism","09164 Neurodegenerative disease",
                                "09191 Unclassified: metabolism","09105 Amino acid metabolism",
                                "09102 Energy metabolism","09131 Membrane transport",
                                "09172 Infectious disease: viral","09143 Cell growth and death",
                                "09171 Infectious disease: bacterial","09108 Metabolism of cofactors and vitamins",
                                "09141 Transport and catabolism","09194 Poorly characterized",
                                "09152 Endocrine system","09145 Cellular community - prokaryotes"),
                      col = c("#FAFD7CFF","#82491EFF","#24325FFF","#B7E4F9FF","#FB6467FF","#526E2DFF","#E762D7FF",
                              "#E89242FF","#FAE48BFF","#A6EEE6FF","#917C5DFF","#00468BFF","#ED0000FF","#42B540FF",
                              "#0099B4FF","#925E9FFF","#FDAF91FF","#AD002AFF")
                        
  ) %>%
    dplyr::slice(1:paln)
  
  ranks <- paste0("rank",c(1:rn))
  rlast <- ranks[length(ranks)]
  
  if (!all(ranks %in% names(keggtax))) {
    stop("Error: need to have taxon levels: ",paste(ranks, collapse = ", "))
  }
  
  kdict <- keggtax[, ranks] %>% distinct()
  
  # start out as all gray
  
  others <- kdict %>%
    filter(! rank2 %in% basecols$rank2 ) 
  
  greyN <- others %>%
    nrow()
  
  greycols <- shades("gray", variation = 0.25, ncolor = greyN)
  othersc <- others %>%
    mutate(greycol = greycols)
  
  kpal <- kdict %>%
    left_join(basecols, by = c("rank2")) %>%
    group_by(rank2, col) %>%
    summarise(catN = n(),
              !!rlast := list(.data[[rlast]])) %>%
    mutate(cols = ifelse(is.na(col),
                         list(),
                         list(shades(col, ncolor = catN, variation = 0.35)))) %>%
    unnest(c("cols",rlast)) %>%
    left_join(othersc) %>%
    mutate(cols = if_else(is.na(cols), greycol, cols)) %>%
    select(rank2,all_of(rlast),cols)
  
  return(structure(kpal$cols,names = kpal[[rlast]]))
  
}
