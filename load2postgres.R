setwd("~/OneDrive - The University of Chicago/shotgun/20230629_SSG76/")

library(RPostgreSQL)
library(tidyverse)

killDbConnections <- function (){

  all_cons <- dbListConnections(dbDriver("PostgreSQL"),
                                host="128.135.41.72",dbname="clinical_db",
                                user="dfi_admin",password="dfibugs")

  print(all_cons)

  for(con in all_cons)
    +  dbDisconnect(con)

  print(paste(length(all_cons), " connections killed."))

}

killDbConnections()

con <- dbConnect(dbDriver("PostgreSQL"),
                 host="128.135.41.72",dbname="clinical_db",
                 user="dfi_admin",password="dfiadmin2022")
dbListTables(con)

# trim stats --------------------------------------------------------------
stats <- readRDS("trimStats.rds")

orgstats <- stats %>%
  separate(seq_id, c("dir","foodir","seq_id"), sep = "/") %>%
  select(-dir, -foodir) %>%
  filter(grepl("[UCMIULDHTF]{2,}",seq_id),
         !grepl("Zhenrun", seq_id))

# load all trimstats: 21
orgstats %>%
  dplyr::count(seq_id) %>%
  view

dbWriteTable(con, "trimstats", orgstats, row.names = FALSE, append = TRUE)

# ID lookup ---------------------------------------------------------------
lkidcheck <- orgstats %>%
  distinct(seq_id) %>%
  mutate(ID = gsub("^[0-9]{6}-","", seq_id),
         ID = if_else(grepl("[HTCUMIFLD]+-[0-9]{3}[-.][0-9]{2}", ID),
                      ID,
                      NA_character_)) %>%
  separate(ID, c("study","subj","tp")) %>%
  mutate(ID = paste0(study, "-", subj, ".", tp),
         ID = if_else(is.na(study),
                      gsub("^[0-9]+-","",gsub("\\.","00",gsub("[0-9]{6}\\-Donor\\-|\\.S","", seq_id))),
                      ID)) %>%
  dplyr::select(seq_id, ID)

# lkidcheck %>% dplyr::slice(c(57, 98))
lkidcheck %>% filter(is.na(ID))

lk <- lkidcheck 

lk

# peek
tbl(con, "lookup") %>%
  collect() %>%
  tail()

dbWriteTable(con, "lookup", lk, row.names = FALSE, append = TRUE)

# kraken2 -----------------------------------------------------------------

kl <-  readRDS("all.kraken2.standard.rds")

clkl <- kl %>%
  as_tibble() %>%
  select(seq_id, uc, taxon, taxid, numseqs, filename) %>%
  mutate(seq_id = gsub("_R1","", seq_id)) %>%
  filter(seq_id %in% lk$seq_id)

clkl %>%
  dplyr::count(seq_id) %>% view

dbWriteTable(con, "kraken2readlevel", clkl, row.names = FALSE, append = TRUE)

# bracken -----------------------------------------------------------------

bra <-  lapply(list.files(pattern = "all.bracken.rds"), readRDS) %>%
  bind_rows() %>%
  dplyr::rename(taxid = taxonomy_id) %>%
  relocate(seq_id) %>%
  mutate(seq_id = gsub("_R1$","",seq_id))

subbra <- bra %>%
  filter(seq_id %in% lk$seq_id) %>%
  as_tibble()

subbra %>%
  dplyr::count(seq_id)

subbra %>%
  dplyr::count(taxonomy_lvl)

dbWriteTable(con, "bracken", subbra, row.names = F, overwrite = F, append = T)

# metaphlan:  ---------------------------------------------------------------

mpa <-  lapply(list.files(pattern = "all.metaphlan.rds"), readRDS) %>%
  bind_rows() %>%
  relocate(seq_id) %>%
  mutate(seq_id = gsub("_R1$","",seq_id),
         taxid = as.numeric(taxid))

submpa <- mpa %>%
  filter(seq_id %in% lk$seq_id) %>%
  as_tibble()

submpa %>%
  dplyr::count(seq_id)

# dbWriteTable(con, "metaphlan", submpa, row.names = F, overwrite = T, append = F)
dbWriteTable(con, "metaphlan", submpa, row.names = F, 
             overwrite = F, append = T)

# tbna = "metaphlan"
#
# dbSendStatement(con, paste0("GRANT SELECT ON ",tbna," TO dfi_lab"))
# dbSendStatement(con, paste0("GRANT SELECT ON ",tbna," TO dfi_user"))
# dbSendStatement(con, paste0("GRANT SELECT ON ",tbna," TO dfi_biobank"))


# targeted ----------------------------------------------

tin <- lapply(list.files(pattern = "targeted.dnmd.rds"), readRDS) %>%
  bind_rows()

tin %>%
  select(-dbsource) %>%
  # left_join(tardict) %>%
  dplyr::count(seq_id, wt = mappedReads, sort = T) %>%
  arrange(seq_id)

tcl <- tin %>%
  # select(-dbsource) %>%
  # left_join(tardict) %>%
  mutate(seq_id = gsub("_R1","", seq_id)) %>%
  left_join(orgstats %>%
              filter(grepl("final pair",countType)) %>%
              select(seq_id, total = count)) %>%
  filter(seq_id %in% lk$seq_id) %>%
  mutate(cpm = mappedReads/total * 1e6,
         rpkm = mappedReads/((slen/1000) * total) * 1e6) %>%
  select(seq_id,sseqid, mappedReads, cpm, rpkm, slen)

tcl %>%
  dplyr::count(seq_id) %>%
  view

dbWriteTable(con, "targeted", tcl, row.names = FALSE, append = TRUE)

# emapper: -------------------------------------------------------

dbListTables(con)

tbl(con,"emapper")

em <- lapply(list.files(pattern = "all.emapper.rds"), readRDS) %>%
  bind_rows()

em %>%
  dplyr::count(sample)

akc <- em %>%
  mutate(seq_id = gsub("_R1\\.$","",sample)) %>%
  # count(seq_id) %>%
  filter(seq_id %in% lk$seq_id) %>%
  select(-sample) %>%
  select(seq_id, everything()) %>%
  mutate(ko_id = str_split(KEGG_ko,",")) %>%
  unnest(ko_id)

akc <- akc %>%
  select(seq_id, ko_id, query_name, starts_with("KEGG"),everything())

# checks
akc %>%
  dplyr::count(seq_id) %>%
  view

dbWriteTable(con, "emapper", akc, row.names = FALSE, append = TRUE)

dbDisconnect(con)
