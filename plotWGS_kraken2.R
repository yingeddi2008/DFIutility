argoptions <- commandArgs(trailingOnly = TRUE)
runid <- argoptions[1]

library(RPostgreSQL)
library(tidyverse)

con <- dbConnect(dbDriver("PostgreSQL"),
                 host = "128.135.41.183",
                 dbname="dfi_commensal_library",
                 use="ericlittmann",
                 password="dfibugs")

dbListTables(con)

tax <- tbl(con, "kraken2_contigs") %>%
  filter(grepl(runid,seq_id)) %>%
  collect() 

taxmap <- tbl(con, "ncbi_blast_taxdump") %>%
  collect()

chgUnclassified <- function(x) {if_else(x == "", "unclassified", x)}

t <- tax %>%
  mutate(taxfoo = str_extract(taxon, "taxid [0-9]+")) %>%
  separate(taxfoo, c("taxidfoo","tax_id"), convert = T) %>%
  left_join(taxmap) %>%
  add_count(seq_id, wt = length, name = "totalLen") %>%
  mutate(pctseqs = length/totalLen) %>%
  count(seq_id, kingdom, phylum, class, order, family, genus, species, 
        wt = pctseqs,
        name = "pctseqs") %>%
  rename_at(vars(kingdom:species), tools::toTitleCase) %>%
  mutate_at(vars(Kingdom:Species), chgUnclassified) %>%
  mutate(genLab = Genus,
         Genus = paste0(Phylum, "-", Class, "-", Order, "-", Family, "-", Genus))

source("~/OneDrive - The University of Chicago/DFIutility/getRdpPal.R")

taxpal <- getRdpPal(t)

scales::show_col(unique(taxpal))  

t %>%
  arrange(Kingdom, Phylum, Class, Order, Family, Genus) %>% 
  mutate(Genus = factor(Genus, levels = unique(Genus))) %>% 
  group_by(seq_id) %>% 
  arrange(Genus) %>% 
  mutate(cum.pct = cumsum(pctseqs), 
         y.text = (cum.pct + c(0, cum.pct[-length(cum.pct)]))/2) %>% 
  ungroup() %>%
  dplyr::select(-cum.pct) %>% 
  mutate(#tax.label=ifelse(Species=="unclassified",paste(Family,Genus,Species,sep="\n"),paste(Genus,Species,sep="\n")),
    tax.label=ifelse(grepl("unclassified$",Genus),genLab,gsub(" ","\n",Species)), 
    tax.label = ifelse(pctseqs >= .1, as.character(tax.label), "")) %>%
  ggplot(aes(seq_id, pctseqs)) +
  theme_bw() +
  geom_col(aes(fill = Genus)) +
  scale_fill_manual(values = taxpal) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_text(aes(y=1-y.text,x=seq_id,label=tax.label),lineheight=0.6,size=2.5, angle = 90) +
  labs(y = "contig length % abundance", title = "kraken2 taxonomy") 

ggsave(paste0(runid,".isolate.kraken2.pdf"), height = 6, width = 8.50)
