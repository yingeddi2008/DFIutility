# generate nanopore long and Illumina reads length distribution and hybrid assembly results

library(stringr)
library(tidyverse)
library(lubridate)

argoptions <- commandArgs(trailingOnly = TRUE)
nanopath <- argoptions[1]
nanodir <- basename(nanopath)

if (is.na(argoptions[2])){
  outdir <- getwd()
} else {
  outdir <- argoptions[2]
}


calNx <- function(lenArr, nx = 0.5){
  len.sorted <- rev(sort(lenArr))
  len.sorted[cumsum(len.sorted) >= sum(len.sorted)*nx][1]
}
calWidth <- function(nsamples){ return(5 + 1.3*nsamples) }

# nanopath <- "/Volumes/pamer-lab/DFI_MMF/nanopore/Bioinformatics/20200922"

# all nanopore long reads
nanofais <- Sys.glob(file.path(nanopath,"*/*.fastq.fai"))
nanofais <- nanofais[grepl("barcode[0-9]{2}.fastq.fai",nanofais)]

nanotmplen <- NULL

for (f in nanofais){
  
  iso <- str_split(f,"/")[[1]][8]
  
  nanotmplen[[iso]] <- read_tsv(f, col_names = F) %>%
    select(X2) %>%
    mutate(isolate = iso)
  
}

nanostats <- bind_rows(nanotmplen) %>%
  group_by(isolate) %>%
  summarise(median = median(X2),
            mean = mean(X2),
            N50 = calNx(X2),
            max = max(X2)) %>%
  mutate(seqType = "Nanopore")

# Illumina 
illufais <- Sys.glob(file.path(nanopath,"*/*/*.fai"))
illufais <- illufais[grepl("hybridAssemblyTrim",illufais)]

illutmplen <- NULL

for (f in illufais){
  
  iso <- str_split(f,"/")[[1]][8]
  
  illutmplen[[iso]] <- read_tsv(f, col_names = F) %>%
    select(X2) %>%
    mutate(isolate = iso)
  
}

illustats <- bind_rows(illutmplen) %>%
  group_by(isolate) %>%
  summarise(median = median(X2),
            mean = mean(X2),
            N50 = calNx(X2),
            max = max(X2)) %>%
  mutate(seqType = "Illumina")

# summary
sumfn <- Sys.glob(file.path(nanopath,"summary*txt"))

summ <- read_table(sumfn,
                   col_names = F) %>%
  mutate(path=if_else(grepl("^##",X1), X1, NA_character_)) %>%
  fill(path) %>%
  filter(!grepl(X1,pattern="^##")) %>% 
  separate(X1, c("contigNo","length","depth","circular"), sep = " ") %>%
  mutate(contigNo = as.numeric(gsub(">","",contigNo)),
         length = as.numeric( gsub("length=","",length) ),
         depth = as.numeric( gsub("depth=|x","",depth) ),
         circular = if_else(grepl("true",circular), "Yes","No"),
         isolate = sapply(str_split(path,"\\/"), function(x) x[2]),
         guppy = if_else(grepl("guppy_fastq", path), "Yes", "No"),
         type = if_else(grepl("hybridAssemblyTrim",path),"moretrimmed","noTrim")) 

goodIsos <- summ %>%
  group_by(path) %>%
  summarise(guppy = unique(guppy),
            isolate = unique(isolate),
            type = unique(type),
            nContigs = n(),
            totalLength = sum(length),
            nCircular = sum(circular == "Yes"),
            goodAssembly = case_when(
              nCircular >= 1 & nContigs <= 5 ~ "Yes",
              nCircular >= 1 & (nContigs - nCircular) <= 2  & nContigs < 8 ~ "Yes",
              TRUE ~ "No")) %>%
  filter(goodAssembly == "Yes", !is.na(isolate)) %>%
  distinct(isolate) %>%
  `$`(isolate)

# combine both and plot ---------------------------------------------------

stats <- bind_rows(nanostats, illustats)

bind_rows(nanotmplen) %>%
  mutate(seqType = "Nanopore") %>%
  bind_rows(  bind_rows(illutmplen) %>%
                mutate(seqType = "Illumina") ) %>%
  mutate(goodAssembly = if_else(isolate %in% goodIsos, "Good","Nah")) %>%
  ggplot(aes(X2)) +
  geom_histogram(aes(fill = goodAssembly), bins = 120 ) +
  geom_vline(aes(xintercept = median,
                 color = "median"), data = stats,linetype = 2) +
  geom_vline(aes(xintercept = N50,
                 color = "N50"), data = stats, linetype = 2) +
  geom_vline(aes(xintercept = max,
                 color = "max"), data = stats, linetype = 2) +
  facet_grid(seqType ~ isolate, scales = "free") +
  theme_bw() +
  scale_x_log10() +
  # scale_y_continuous(trans = "log2") +
  scale_fill_brewer(palette = "Set1") +
  scale_color_manual(name = "statistics", 
                     values = c(median = "#d17b11", 
                                max = "#19cf19", 
                                N50 = "#d111a7")) +
  labs(x = "read length", y = "# of reads", 
       fill = "hybrid assembly" ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(size = 6.5))

wd <- calWidth(nrow(nanostats))

ggsave(file.path(outdir,
                 paste0("nanopore.isolates.hist.",
                        nanodir,".pdf")),height = 5.3, width = wd)
