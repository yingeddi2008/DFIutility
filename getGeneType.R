getGeneType <- function(seq_ids,prokka_table) {
  
  totIsos <- length(seq_ids)
  
  print("any gene with a product assignment...")
  print(paste("finding core genome for",totIsos,"isolates..."))
  
  tmp <- prokka_table %>%
    filter(seq_id %in% seq_ids,
           !(product %in% c("hypothetical protein","putative protein"))) %>%
    dplyr::count(seq_id, product, name = "isoPdtCnt")  %>%
    add_count(seq_id, wt = isoPdtCnt, name = "totalGenes")

  geneTyped <- tmp %>%
    dplyr::count(product, name = "pdtCnt") %>%
    mutate(gene_type = case_when(
      pdtCnt == totIsos ~ "Core",
      pdtCnt == 1 ~ "Unique",
      TRUE ~ "Shared"
    ))
  
  return(tmp %>% left_join(geneTyped))
  
}
