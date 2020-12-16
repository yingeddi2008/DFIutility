# read in prokka table for WGS runs

source("~/Documents/Eddi/DFIutility/readin_prokka_gff.R")
source("~/Documents/Eddi/DFIutility/readin_prokka_tables.R")
source("~/Documents/Eddi/DFIutility/readin_kraken2_contigs.R")

library(RPostgreSQL)
library(tidyverse)
library(data.table)

argopts <- commandArgs(trailingOnly = TRUE) 
path = argopts[1]
con <- dbConnect(dbDriver("PostgreSQL"),
                 host = "128.135.41.183",
                 dbname="dfi_commensal_library",
                 use="ericlittmann",
                 password="dfibugs")

# dbListTables(con) %>% sort()
ldF = as.logical(argopts[2])
# example path
# path = "/mnt/Data/sdc/WGS/201105_M00621_0042_000000000-JD5L6-EP-10S"

# read in prokka annotations and sequences --------------------------------

ptbls <- readin_prokka_tables(file.path(path,"prokka"))

# load annotations --------------------------------------------------------
# tbl(con, "prokka_annotations")

pann <- ptbls[[1]] %>%
  tibble() %>%
  filter(ftype != "gene") %>%
  separate(filename, c("seq_id","trashfn"), remove = F, sep = "/") %>%
  mutate(trashfn = NULL) 

if (nrow(pann) > 0 & ldF){
  dbWriteTable(con, "prokka_annotations", pann, row.names = F, append = T)
} else {
  warning(">> Nothing to write for annotation!!")
}

# load sequences ----------------------------------------------------------
# tbl(con, "prokka_sequences")

pseqs <- ptbls[[2]] %>%
  tibble() %>%
  select(-filename)

if (nrow(pseqs) > 0 & ldF){
 dbWriteTable(con, "prokka_sequences", pseqs, row.names = F, append = T)
} else {

 warning(">> Nothing to write for sequences!!")
}
# read in gff -------------------------------------------------------------

gffcols <- tbl(con, "prokka_annotations_location") %>%
  colnames()

pgff <- readin_prokka_gff(file.path(path,"prokka"))

cgff <- pgff %>%
  tibble() %>%
  filter(type != "gene") %>%
  separate(filename, c("seq_id","trashfn"), remove = F, sep = "/") %>%
  # dplyr::count(type)
  select(all_of(gffcols))

if (nrow(cgff) > 0 & ldF){
  dbWriteTable(con, "prokka_annotations_location", cgff, row.names = F, append = T)
} else {
  warning(">> Nothing to write for gff information!!")
}


# look up table -----------------------------------------------------------
lkcls <- tbl(con, "lookup_table") %>%
  colnames()

morelk <- tibble(seq_id = unique(pann$seq_id)) %>%
  mutate(interid = gsub("^EP-[0-9]+[a-zA-Z]?-","", seq_id),
         msk_id = gsub("-|-0",".",interid),
         interid = gsub("-[0-9]+$","",interid)) %>%
  separate(interid, c("pre","num"), convert = T) %>%
  mutate(num = formatC(num, width = 4, flag = "0"),
         donor_id = paste0(pre,num),
         pre = NULL, num = NULL,
         cmic_id = NA_character_,
         genome_id = NA_real_) %>%
  select(all_of(lkcls))

morelk

# if (nrow(morelk) > 0 & ldF ){
#  dbWriteTable(con, "lookup_table", morelk, row.names = F, append = T)
#} else {
#  warning(">> Nothing to write for lookup table!!")
#}

# every 16S gene annotation ---------------------------------------------
rnaref <- "/home/pamerlab/Documents/Eddi/refseq_rna/refseq_rna"

seqs16s <- pann %>%
  filter(grepl("16S riboso", product)) %>%
  left_join(pseqs) 

map16s <- seqs16s %>%
  select(seq_id, locus_tag)

seqs16s %>%
  transmute(fasta = paste0(seq_name,"\n",
                           nuc_sequence)) %>%
  write.table("temp.16S.fasta", col.names = F, row.names = F, quote = F)

tax <- tbl(con, "ncbi_blast_taxdump") %>% collect()

system(paste0("/home/pamerlab/Downloads/ncbi-blast-2.10.0+/bin/blastn -query temp.16S.fasta -db ",
              rnaref, " -outfmt '6 std score qlen qcovs sstrand slen staxid' -max_target_seqs 30 -num_threads 8",
              " -out temp.16S.blastn"))

blast <- readr::read_delim("temp.16S.blastn", delim = "\t", 
                           col_names = c("qaccver", "saccver", "pident", "length", 
                                         "mismatch", "gapopen", "qstart", "qend", 
                                         "sstart", "send","evalue", "bitscore", "score", 
                                         "qlen", "qcovus", "sstrand", "slen","staxid")) %>%
  mutate(qcov = (qend - qstart + 1)/qlen * 100,
         iterPos = ifelse(sstrand == "minus", send, sstart),   #switch start and end position in minus
         send = ifelse(sstrand == "minus", sstart, send),
         sstart = ifelse(sstrand == "minus", iterPos, sstart),
         iterPos = NULL) %>%
  group_by(qaccver) %>%
  top_n(5, wt = bitscore) %>%
  left_join(tax, by = c("staxid" = "tax_id")) %>%
  left_join(map16s, by = c("qaccver" = "locus_tag")) %>%
  write_csv("16S.top5.taxonomy.csv")
  
# load kraken2 -------------------------------------------
kraken2 <- readin_kraken2_contigs(file.path(path,"kraken2"))
if(nrow(kraken2) >0 & ldF){
  print("appending kraken2 contigs")
  dbWriteTable(con,"kraken2_contigs",kraken2,row.names=F,append=T)
}else{
  print("no kraken2 contigs to append")
}

# disconnect postgres ------------------------------------
if (dbDisconnect(con)){
  warning(">> Disconnected from PostgreSQL database!")
}

# system("rm temp.16S.fasta temp.16S.blastn")
