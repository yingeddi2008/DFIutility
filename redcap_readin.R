rm(list=ls())
library(tidyverse)

redcap_lt <- function(directory){
  # directory <- "~/Box/Erics_Stuff/Clinical_Data/20210709_redcap/"
  file_events <-
    dir(path = directory,
        pattern = "MicrobiomeInLiverTra-EventTableBioinforma_DATA",
        full.names = T)
  file_samples <-
    dir(path = directory,
        pattern = "MicrobiomeInLiverTra-StoolSampleIDsCollec_DATA",
        full.names = T)
  
  lt_samples <- read.csv(file = file_samples,
                         stringsAsFactors = F) %>%
    mutate(studyid = gsub(",.+", "", record_id)) %>%
    arrange(studyid, mrn) %>%
    fill(mrn) %>%
    select(studyid,
           mrn,
           sampleid = sample_id,
           date_collected,
           sampletype = sample_type) %>%
    filter(sampleid != "") %>%
    mutate(db = "LiverTransplant",
           date_collected = as.Date(date_collected))
  
  lt_events <- read.csv(file = file_events,
                        stringsAsFactors = F) %>%
    fill(mrn) %>%
    filter(event_date != "",
           grepl(pattern = "Transplant", x = event_description)) %>%
    mutate(studyid = gsub(",.+", "", record_id)) %>%
    select(studyid, mrn, date_transplant = event_date) %>%
    unique() %>%
    mutate(db = "LiverTransplant",
           date_transplant = as.Date(date_transplant)) %>% 
    distinct()
  
  lt_int <- full_join(lt_samples,lt_events)
  
  return(lt_int)
}

redcap_micu <- function(directory){
  # directory <- "~/Box/Erics_Stuff/Clinical_Data/20210709_redcap/"
  file_micu <-
    dir(path = directory,
        pattern = "MICUMicrobiomeDataba-MRNStoolSampleCollec_DATA",
        full.names = T)
  
  micu <- read.csv(file = file_micu,
                   stringsAsFactors = F) %>%
    dplyr::rename(date_collected = stool) %>%
    mutate(studyid = ifelse(studyid == "", NA, studyid)) %>%
    fill(mrn, studyid) %>%
    filter(date_collected != "") %>%
    mutate(date_collected = (as.Date(date_collected))) %>%
    select(studyid, mrn, sampleid, date_collected) %>%
    mutate(db = "MICU") %>%
    mutate(sampletype = "Stool")
  
  return(micu)
}

redcap_ht <- function(directory){
  # directory <- "~/Box/Erics_Stuff/Clinical_Data/20210709_redcap/"
  file_ht <-
    dir(path = directory,
        pattern = "MicrobiomeDSA-NewSampleIDsTranspla_DATA",
        full.names = TRUE)
  
  ht <- read.csv(file = file_ht,
                 stringsAsFactors = F) %>%
    dplyr::rename(studyid = study_id,
                  date_transplant = ohtdate,
                  sampleid = sample_id) %>%
    mutate(
      studyid = ifelse(studyid == "", NA, studyid),
      date_transplant = ifelse(date_transplant == "", NA, date_transplant),
      date_transplant = ifelse(
        !is.na(mrn) & is.na(date_transplant),
        "no date",
        date_transplant
      ),
      date_transplant = (as.Date(date_transplant))
    ) %>%
    fill(mrn, studyid, date_transplant) %>%
    filter(!sampleid == "") %>%
    select(studyid, mrn, sampleid, date_collected, date_transplant) %>%
    mutate(
      db = "HeartTransplant",
      date_collected = as.Date(date_collected),
      sampletype = "Stool"
    )
  
  return(ht)
}

redcap_ld <- function(directory){
  # directory <- "~/Box/Erics_Stuff/Clinical_Data/20210709_redcap/"
  file_ld <-
    dir(path = directory,
        pattern = "LiverDiseaseMicrobio-SampleIDMRNBioinform_DATA",
        full.names = TRUE)
  
  ld <- read.csv(file = file_ld,
                 stringsAsFactors = F) %>% 
    dplyr::rename(studyid = subject_id,
                  date_collected = collection_date,
                  sampleid = sample_id,
                  sampletype = sample_type) %>%
    mutate(
      studyid = ifelse(studyid == "", NA, studyid),
      date_collected = ifelse(date_collected == "", NA, date_collected),
      # date_transplant = ifelse(
      #   !is.na(mrn) & is.na(date_transplant),
      #   "no date",
      #   date_transplant
      # ),
      date_collected = (as.Date(date_collected))
    ) %>%
    fill(mrn, studyid) %>%
    filter(!sampleid == "") %>%
    select(studyid,sampletype, mrn, sampleid, date_collected) %>%
    mutate(
      db = "LiverDisease",
      date_collected = as.Date(date_collected),
      sampletype = ifelse(sampletype == 1, "Stool", "Other")
    )
  
  return(ld)
}

redcap_wrapper <- function(directory){
  lt <- redcap_lt(directory)
  micu <- redcap_micu(directory)
  ht <- redcap_ht(directory)
  ld <- redcap_ld(directory)
  total <- lt %>% 
    bind_rows(micu, ht, ld) %>% 
    mutate(sampleid = gsub(pattern = " ", replacement = "", sampleid)) %>% 
    separate(sampleid, c("LD", "subn", "tp"), convert = T, remove = F) %>%
    mutate(ID = paste0(LD, "-", formatC(subn, width = 3, flag = "0"),".",
                       formatC(tp, width = 2, flag = "0"))) %>% 
    select(-c(LD,subn,tp))
  
  return(total)
}

redcap_tables <- redcap_wrapper(directory = "~/Box/Erics_Stuff/Clinical_Data/20210709_redcap/")

# Establish connection with Postgres
con <- dbConnect(dbDriver("PostgreSQL"),
                 host="128.135.41.32",
                 dbname="clinical_db",
                 user="dfi_admin",
                 password="dfibugs")

# Load most recent postgres lookup table (to shotgun sequences)
lookup <- tbl(con,"lookup") %>% collect()

bile <- tbl(con, "bile_v3") %>% collect()

scfa <- tbl(con, "scfa_v3") %>% collect()

shotgun_lookup <- redcap_tables %>% 
  full_join(lookup, by = "ID") %>% 
  dplyr::rename(shotgunSeq_id = seq_id)

# load most recent 16s phyloseq object
phy <- readRDS("/Volumes/pamer-lab/DFI_MMF/MiSeq/Bioinformatics/dada2rds/finalPhy_rdp.merged.wMeta.rds")

subphy <- subset_samples(phy, grepl("_[MICU\\.\\_\\-]{3,}[0-9]{2}", sample_names(phy)) | 
                           grepl("_[UC\\.\\_\\-]{3,}[0-9]{2}",  sample_names(phy)) | 
                           grepl("_[HT\\.\\_\\-]{3,}[0-9]{2}",  sample_names(phy)) |
                           grepl("_[LD\\.\\_\\-]{3,}[0-9]{2}", sample_names(phy)))

# tmp16s <- as_tibble(sample_names(subphy))

subphy2 <- subset_samples(subphy, !grepl("CC",sample_names(subphy)))

phy_lookup <- sample_data(subphy2) %>%
  data.frame() %>%
  rownames_to_column(var = "MiSeq_id") %>%
  as_tibble() %>%
  select(MiSeq_id, reads.in, reads.out) %>%
  separate(MiSeq_id, c("pool","ID"), sep = "_", remove = F) %>%
  group_by(ID) %>%
  top_n(1, wt = reads.out) %>%
  separate(ID, c("study","subj","tp")) %>%
  mutate(ID = paste0(study, "-",subj, ".", tp)) %>%
  select(MiSeq_id, ID)

# join shotgun data with 16s data
redcap_genomics <- phy_lookup %>% 
  full_join(shotgun_lookup) 

# join metabolomics (bile acid and pfbbr) and clean-up some errant metabolomicsIDs
metab <- tibble(metabolomicsID = unique(c(bile$metabolomicsID,scfa$metabolomicsID))) %>% 
  mutate(metabolomicsID_corrected = gsub(pattern = " ", replacement = "", metabolomicsID),
         metabolomicsID_corrected = gsub(pattern = "[a-z]+$", replacement = "", metabolomicsID_corrected),
         metabolomicsID_corrected = case_when(metabolomicsID_corrected == "HT_005_01_RI"~"HT_005_01",
                                              metabolomicsID_corrected == "U_C_003_05"~"UC_003_05",
                                              metabolomicsID_corrected == "U_C_003_06"~"UC_003_06",
                                              metabolomicsID_corrected == "U_C_003_07"~"UC_003_07",
                                              metabolomicsID_corrected == "U_C_006_01"~"UC_006_01",
                                              metabolomicsID_corrected == "U_C_005_01"~"UC_005_01",
                                              metabolomicsID_corrected == "U_C_007_02"~"UC_007_02",
                                              metabolomicsID_corrected == "U_C_004_03"~"UC_004_03",
                                              metabolomicsID_corrected == "MICU_048_01_RI"~"MICU_048_01",
                                              TRUE ~ metabolomicsID_corrected)
         ) %>% 
  separate(metabolomicsID_corrected, c("LD", "subn", "tp"), remove = F, convert = T) %>%
  mutate(ID = paste0(LD, "-", formatC(subn, width = 3, flag = "0"),".",
                     formatC(tp, width = 2, flag = "0"))) %>% 
  select(metabolomicsID, ID)

# join all metabolomics and genomics lookups and clean up some syntax
redcap_genomics_metab <- redcap_genomics %>%
  full_join(metab, by = "ID") %>% 
  mutate(sampletype = tolower(sampletype))

#  Re-write table to redcap
dbWriteTable(con, "redcap_tbl_v3", redcap_genomics_metab,
             row.names = F, append = F, overwrite = T)
dbSendStatement(con, "GRANT SELECT ON redcap_tbl_v3 TO dfi_lab")
dbSendStatement(con, "GRANT SELECT ON redcap_tbl_v3 TO dfi_user")

dbDisconnect(con)

# load redcap table
# temp <- tbl(con, "redcap_tbl_v3") %>% collect()



