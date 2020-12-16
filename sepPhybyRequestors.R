library(tidyverse)
library(dplyr)
library(phyloseq)

calWidth <- function(nsamples){ return(5 + 0.22*nsamples) }

argopts <- commandArgs(trailingOnly = TRUE)
inFilePath <- argopts[1]
# inFilePath <- "/Volumes/pamer-lab/DFI_MMF/MiSeq/Bioinformatics/runPools/200921_M00621_0028_000000000-D9FVY-EP-SM32/SM32"
poolId <- basename(inFilePath)

if (is.na(argopts[2])){
  blF <- FALSE
} else {
  blF <- as.logical(argopts[2])
}

# rdp plot by requester -----------------------------------------------------
phy <- readRDS(file.path(inFilePath,"cleanFastq/finalPhy_rdp.rds"))

# blast tax read in -------------------------------------------------------
if (blF){
  bln <- read_csv(file.path(inFilePath, "cleanFastq","asv_seqs.blastn.tax"))
  
  btax <- bln %>%
    group_by(qaccver) %>%
    top_n(1, wt = bitscore) %>%
    top_n(1, wt = pident) %>%
    top_n(1, wt = (qend - qstart)/qlen) %>%
    top_n(-1, wt = evalue) %>% 
    add_count(species, name = "spwt") %>%
    sample_n(1, weight = spwt) %>% 
    dplyr::rename(Superkingdom = kingdom,
                  Species = species,
                  Genus = genus,
                  Family = family,
                  Order = order,
                  Class = class,
                  Phylum = phylum) %>%
    ungroup()
  
  phytax <- tax_table(as.matrix(btax[,2:ncol(btax)]))
  taxa_names(phytax) <- btax$qaccver
  
  blastphy <- phy
  tax_table(blastphy) <- phytax
  
  #  saveRDS(blastphy, "finalPhy_blast.rds")
  saveRDS(blastphy, file.path(inFilePath,"cleanFastq","finalPhy_blast.rds"))  
}


# separate by requestor ---------------------------------------------------

reqtors <- unique(sample_data(phy)$request.by)

for (rh in reqtors){
  # temp phyloseq object
  tempphy <- subset_samples(phy, request.by == rh)
  tempphy <- subset_taxa(tempphy, taxa_sums(tempphy) > 0)
  # write out a csv format
  tax_table(tempphy) %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column(var = "taxonID") %>%
    write_csv(paste0(poolId,"_",gsub(" ","",rh),".taxTable_rdp.csv"))
  
  if (blF){
    tempblastphy <- subset_samples(blastphy, request.by == rh)
    tempblastphy <- subset_taxa(tempblastphy, taxa_sums(tempblastphy) > 0)
  
  # write out a csv format
  tax_table(tempblastphy) %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column(var = "taxonID") %>%
    write_csv(paste0(poolId,"_",gsub(" ","",rh),".taxTable_blast.csv"))
  
  saveRDS(tempblastphy,paste0(poolId,"_", gsub(" ","",rh),".finalPhy_blast.rds"))
  }
  
  otu_table(tempphy) %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column(var = "sampleID") %>%
    write_csv(paste0(poolId,"_",gsub(" ","",rh),".otuTable_rdp.csv"))
  
  sample_data(tempphy) %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column(var = "sampleID") %>%
    write_csv(paste0(poolId,"_",gsub(" ","",rh),".metaTable_rdp.csv"))
  
  Biostrings::writeXStringSet(refseq(tempphy),
                  paste0(poolId,"_", gsub(" ","",rh),".asv_seqs.fasta"))
  
  saveRDS(tempphy,paste0(poolId,"_", gsub(" ","",rh),".finalPhy_rdp.rds"))
  
}
