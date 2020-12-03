library(tidyverse)
library(dplyr)
library(yingtools2)
library(data.table)
library(grid)
library(gridExtra)
library(phyloseq)
library(stringr)

source("/Volumes/pamer-lab/Eric.Littmann/R_reference_files/prokka_functions_EL_REF.R")

opts <- commandArgs(trailingOnly = TRUE)

inFilePath <- opts
poolID <- str_split(tail(str_split(inFilePath,"/")[[1]],1),"\\.")[[1]][1]

calWidth <- function(nsamples){ return(3 + 0.28*nsamples) }

# poolID <- "SM36"
# inFilePath <- "/Volumes/pamer-lab/DFI_MMF/MiSeq/Bioinformatics/runPools/201005_M00621_0033_000000000-D9M7G-EP-SM36/SM36/"
outFilePath <- getwd()

# rdp plot ----------------------------------------------------------------
phy <- readRDS(inFilePath)

treads <- sum(sample_data(phy)$reads.in)
freads <- sum(sample_data(phy)$reads.out)

t <- get.otu.melt(phy) %>%
  replace_na(list(Species="unclassified",
                  Genus="unclassified",
                  Family="unclassified",
                  Order="unclassified",
                  Class="unclassified",
                  Phylum="unclassified")) 

pal <- get.yt.palette2(t)

wd <- calWidth(nsamples(phy))

if (length(unique(t$group)) > 1){
  
  totals <- t %>%
    group_by(group,sample) %>%
    summarize(total=sum(numseqs))
  
  gg <- t %>%
    group_by(sample,group,Superkingdom,Phylum,Class,Order,Family,Genus, Species) %>%
    summarize(pctseqs=sum(pctseqs)) %>%
    ungroup() %>%
    arrange(Superkingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
    mutate(Species = factor(Species, levels = unique(Species))) %>% 
    group_by(group, sample) %>% 
    arrange(Species) %>% 
    mutate(cum.pct = cumsum(pctseqs), 
           y.text = (cum.pct + c(0, cum.pct[-length(cum.pct)]))/2) %>% 
    ungroup() %>%
    dplyr::select(-cum.pct) %>% 
    mutate(#tax.label=ifelse(Species=="unclassified",paste(Family,Genus,Species,sep="\n"),paste(Genus,Species,sep="\n")),
      tax.label=ifelse(Genus=="unclassified",paste(Family,Genus,sep="\n"),Genus), 
      tax.label = ifelse(pctseqs >= .1, as.character(Genus), "")) %>%
    left_join(totals) %>%
    ggplot(aes(x=sample,y=pctseqs)) +
    geom_bar(aes(fill=Species),stat="identity") +
    geom_text(aes(y=1-y.text,x=sample,label=tax.label),angle=90,
              lineheight=0.6,size=2.5) +
    scale_fill_manual(values=pal) +
    ggtitle("blast 16s taxonomy") +
    ylab("16S % Abundance") +
    facet_grid(. ~ group, scales = "free",space = "free") +
    xlab("") +
    theme_bw() +
    theme(#axis.text.x=element_text(angle=90,size=12),
      axis.text.x=element_blank(),
      strip.text.x=element_text(size=15),
      legend.position="none") 
  
  gg2 <- totals %>%
    ggplot(aes(x=sample,y=total)) +
    geom_bar(color = "black", alpha = 0.65, fill = "grey", width=0.75, stat="identity") +
    theme_bw() +
    labs(x = "", y = "Number of reads", 
         caption = paste0("Total reads: ",treads,
                          "; high quality reads: ",signif(freads/treads*100,4),
                          "%; # of samples:", nsamples(phy))) +
    facet_grid(. ~ group, scales = "free",space = "free") +
    geom_hline(yintercept=5000,linetype="dashed", color = "blue") +
    geom_hline(yintercept=1000,linetype="dashed", color = "red") +
    theme(axis.text.x=element_text(angle=90, vjust = 0.5, hjust = 1 )) +
    scale_y_continuous(trans = "log10")
  
  pdf(file = file.path(outFilePath, paste0(poolID,"_blast_barplot_ind.pdf")),
      height = 8, width = wd)
  # pdf(file = "foo.pdf", height = 9, width = 7)
  # ggarrange(gg,gg2, ncol = 1)
  
  gstack <- gg.stack(gg,gg2,heights=c(4,1), newpage = FALSE)
  
  dev.off()
  
} else {
  
  totals <- t %>%
    group_by(sample) %>%
    summarize(total=sum(numseqs))
  
  gg <- t %>%
    group_by(sample,Superkingdom,Phylum,Class,Order,Family,Genus, Species) %>%
    summarize(pctseqs=sum(pctseqs)) %>%
    ungroup() %>%
    arrange(Superkingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
    mutate(Species = factor(Species, levels = unique(Species))) %>% 
    group_by(sample) %>% 
    arrange(Species) %>% 
    mutate(cum.pct = cumsum(pctseqs), 
           y.text = (cum.pct + c(0, cum.pct[-length(cum.pct)]))/2) %>% 
    ungroup() %>%
    dplyr::select(-cum.pct) %>% 
    mutate(#tax.label=ifelse(Species=="unclassified",paste(Family,Genus,Species,sep="\n"),paste(Genus,Species,sep="\n")),
      tax.label=ifelse(Genus=="unclassified",paste(Family,Genus,sep="\n"),Genus), 
      tax.label = ifelse(pctseqs >= .1, as.character(Genus), "")) %>%
    left_join(totals) %>%
    ggplot(aes(x=sample,y=pctseqs)) +
    geom_bar(aes(fill=Species),stat="identity") +
    geom_text(aes(y=1-y.text,x=sample,label=tax.label),angle=90,
              lineheight=0.6,size=2.5) +
    scale_fill_manual(values=pal) +
    ggtitle("blast 16s taxonomy") +
    ylab("16S % Abundance") +
    xlab("") +
    theme_bw() +
    theme(#axis.text.x=element_text(angle=90,size=12),
      axis.text.x=element_blank(),
      strip.text.x=element_text(size=15),
      legend.position="none") 
  
  gg2 <- totals %>%
    ggplot(aes(x=sample,y=total)) +
    geom_bar(color = "black", alpha = 0.65, fill = "grey", width=0.75, stat="identity") +
    theme_bw() +
    labs(x = "", y = "Number of reads", 
         caption = paste0("Total reads: ",treads,
                          "; high quality reads: ",signif(freads/treads*100,4),
                          "%; # of samples:", nsamples(phy))) +
    geom_hline(yintercept=5000,linetype="dashed", color = "blue") +
    geom_hline(yintercept=1000,linetype="dashed", color = "red") +
    theme(axis.text.x=element_text(angle=90, vjust = 0.5, hjust = 1 )) +
    scale_y_continuous(trans = "log10")
  
  pdf(file = file.path(outFilePath,paste0(poolID,"_blast_barplot_ind.pdf")),
      height = 8, width = wd)
  # pdf(file = "foo.pdf", height = 9, width = 7)
  # ggarrange(gg,gg2, ncol = 1)
  
  gstack <- gg.stack(gg,gg2,heights=c(4,1), newpage = FALSE)
  
  dev.off()
}


