getKeggPal <- function(keggtax, rn = 2, legend = F) {
  
  require(yingtools2)
  require(tidyverse)
  require(assertive.types)
  require(assertive.base)

# hard code color and rank1 values ----------------------------------------
  
  basecols <- tibble(rank1 = c('09100 Metabolism',
                               '09190 Not Included in Pathway or Brite',
                               '09120 Genetic Information Processing',
                               '09130 Environmental Information Processing',
                               '09140 Cellular Processes',
                               '09160 Human Diseases',
                               '09150 Organismal Systems'),
                     col = c('#FF7F0EFF','#dae0eb','#2CA02CFF',
                             '#9467BDFF','#8C564BFF','#E377C2FF','#17BECFFF'))
  
  baseplt <- basecols %>%
    group_by(rank1) %>%
    mutate(cols = list(shades(col, variation = 0.25, ncolor = 3)),
           order = list(c(1,2,3))) %>%
    unnest(c("cols","order")) 
  
  basepal <- structure(baseplt$cols, names = baseplt$cols)
  
  lg <- baseplt %>%
    ggplot(aes(order, rank1)) +
    geom_tile(aes(fill = cols)) +
    scale_fill_manual(values = basepal) +
    theme_minimal() +
    labs(x = "", y = "") +
    theme(axis.text.x = element_blank(),
          legend.position = "none",
          axis.text.y = element_text(size = 17),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    scale_y_discrete(position = "right")

# legend output ------------------------------------------------------------
  legend <- assert_is_logical(legend)
  if (legend) { 
    return(lg)
  }

# argument checks ---------------------------------------------------------

  rn <- assert_is_a_number(rn)
  
  if ( between(rn, 2,4) ){
    warning(
      switch(as.character(rn),
      "2" = "Good choice of rank level at 2, and it is the default!",
      "3" = "Rank 3 may work! Possibly a bit busy.",
      "4" = "Too many colors! Don't recommend it!"
      )
    )
  } else {
    stop("rn argument is out of range! Please correct! It should be between 2, and 4.")
  }
  
  # switch to rank2
  # basecols <- tibble( rank2 = c("09182 Protein families: genetic information processing",
  #                               "09183 Protein families: signaling and cellular processes",
  #                               "09181 Protein families: metabolism","09132 Signal transduction",
  #                               "09101 Carbohydrate metabolism","09164 Neurodegenerative disease",
  #                               "09191 Unclassified: metabolism","09105 Amino acid metabolism",
  #                               "09102 Energy metabolism","09131 Membrane transport",
  #                               "09172 Infectious disease: viral","09143 Cell growth and death",
  #                               "09171 Infectious disease: bacterial","09108 Metabolism of cofactors and vitamins",
  #                               "09141 Transport and catabolism","09194 Poorly characterized",
  #                               "09152 Endocrine system","09145 Cellular community - prokaryotes"),
  #                     col = c("#FAFD7CFF","#82491EFF","#24325FFF","#B7E4F9FF","#FB6467FF","#526E2DFF","#E762D7FF",
  #                             "#E89242FF","#FAE48BFF","#A6EEE6FF","#917C5DFF","#00468BFF","#ED0000FF","#42B540FF",
  #                             "#0099B4FF","#925E9FFF","#FDAF91FF","#AD002AFF")
  #                       
  # ) %>%
  #   dplyr::slice(1:paln)
  
  ranks <- paste0("rank",c(1:rn))
  rlast <- ranks[length(ranks)]
  
  if (!all(ranks %in% names(keggtax))) {
    stop("Error: need to have taxon levels: ",paste(ranks, collapse = ", "))
  }
  
  kdict <- keggtax[, ranks] %>% distinct()
  
  # start out as all gray
  others <- kdict %>%
    filter(! rank1 %in% basecols$rank1 ) 
  
  greyN <- others %>%
    nrow()
  
  greycols <- shades("gray", variation = 0.25, ncolor = greyN)
  othersc <- others %>%
    mutate(greycol = greycols)
  
  kpal <- kdict %>%
    left_join(basecols, by = c("rank1")) %>%
    group_by(rank1, col) %>%
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
