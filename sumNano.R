# generate nanopore long and Illumina reads length distribution and hybrid assembly results

library(stringr)
library(tidyverse)
library(lubridate)
library(conflicted)

argoptions <- commandArgs(trailingOnly = TRUE)
nanopath <- argoptions[1]
nanopath <- gsub("/$","",nanopath)
nanodir <- basename(nanopath)

if (is.na(argoptions[2])){
  outdir <- getwd()
} else {
  outdir <- argoptions[2]
}

conflict_prefer("n", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")

calNx <- function(lenArr, nx = 0.5){
  len.sorted <- rev(sort(lenArr))
  len.sorted[cumsum(len.sorted) >= sum(len.sorted)*nx][1]
}

calWidth <- function(nsamples){ return(5 + 1.3*nsamples) }

# nanopath <- "/Volumes/pamer-lab/DFI_MMF/nanopore/Bioinformatics/20210920"

# all nanopore long reads ------------------------------------------------
nanofais <- Sys.glob(file.path(nanopath,"*/*.fastq.fai"))
nanofais <- nanofais[grepl("barcode[0-9]{2}.fastq.fai",nanofais)]

nanotmplen <- NULL

for (f in nanofais){

  iso <- str_split(f,"/")[[1]][9]

  nanotmplen[[iso]] <- read_tsv(f, col_names = F) %>%
    select(X2) %>%
    mutate(isolate = iso)

}

nanostats <- bind_rows(nanotmplen) %>%
  group_by(isolate) %>%
  summarise(median = median(X2),
            mean = mean(X2),
            N50 = calNx(X2),
            max = max(X2),
            n=dplyr::n(),
            total_bp = sum(X2)) %>%
  mutate(seqType = "Nanopore")

# Illumina -----------------------------------------------------

illufais <- Sys.glob(file.path(nanopath,"*/*/002*.fai"))
# illufais <- illufais[grepl("hybridAssemblyTrim",illufais)]

illutmplen <- NULL

for (f in illufais){

  iso <- str_split(f,"/")[[1]][9]

  illutmplen[[f]] <- read_tsv(f, col_names = F) %>%
    select(X2) %>%
    mutate(isolate = iso,
           file = f,
           type = if_else(grepl("hybridAssemblyTrim",f),"moretrimmed","noTrim"))

}

illustats <- bind_rows(illutmplen) %>%
  group_by(isolate, file, type) %>%
  summarise(median = median(X2),
            mean = mean(X2),
            N50 = calNx(X2),
            max = max(X2),
            n=dplyr::n(),
            total_bp = sum(X2)) %>%
  mutate(seqType = "Illumina")

# total length of flye and long assembly ----------------------------------

## flye -------------------------------------------------------------------

flye.fns <- Sys.glob(file.path(nanopath,"*/flye/assembly_info.txt"))

flye.tmp <- list()

for (flyfn in flye.fns){
  
  iso <- str_split(flyfn,"/")[[1]][9]
  
  tmpflye <- read_tsv(flyfn) %>%
    transmute(contigNo = `#seq_name`,
           length, 
           circular = factor(circ., levels = c("Y","N"),
                             labels = c("Yes","No")),
           isolate = iso,
           cov = cov.,
           type = "flyeLong") 
  
  denom <- tmpflye %>% 
    slice_max(length, n = 1) %>% 
    pull(cov)
  
  flye.tmp[[iso]] <- tmpflye %>% 
    mutate(depth = cov/denom) %>% 
    select(-cov)
  
}

# bind_rows(flye.tmp) %>% 
#   view

## unicycler long reads only ----------------------------------------

### remove unicycler long reads only assembly

# summary of the assemblies -----------------------------------------
sumfn <- Sys.glob(file.path(nanopath,"summary*txt"))

summ <- read_tsv(last(sumfn),
                   col_names = F) %>%
  mutate(path=if_else(grepl("^##",X1), X1, NA_character_)) %>%
  fill(path) %>%
  filter(!grepl(X1,pattern="^##")) %>%
  separate(X1, c("contigNo","length","depth","circular"), sep = " ") %>%
  mutate(contigNo = gsub(">","",contigNo),
         length = as.numeric( gsub("length=","",length) ),
         depth = as.numeric( gsub("depth=|x","",depth) ),
         circular = if_else(grepl("true",circular), "Yes","No"),
         isolate = sapply(str_split(path,"\\/"), function(x) x[2]),
         guppy = if_else(grepl("guppy_fastq", path), "Yes", "No"),
         type = case_when(
           grepl("hybridAssemblyTrim",path) ~ "moretrimmed",
           grepl("hybridAssembly", path) ~ "noTrim",
           grepl("flye", path) ~ "flyeLong",
           TRUE ~ "long"),
         circular = if_else(grepl("long", type, ignore.case = T),
                            NA_character_,
                            circular)) %>% 
  filter(type != "long")

org.summ <- summ %>%
  filter(type != "flyeLong") %>%
  bind_rows(summ %>%
              filter(type == "flyeLong") %>%
              select(-c(length, circular, depth)) %>%
              left_join(bind_rows(flye.tmp))) %>% 
  arrange(isolate, path, 
          desc(length)) %>% 
  relocate(isolate, type)

write_csv(org.summ %>% 
            select(-guppy), 
          file.path(outdir,paste0("detailed.summary.", nanodir,".nano.csv")))

contigstats <- org.summ %>%
  group_by(path) %>%
  summarise(guppy = unique(guppy),
            isolate = unique(isolate),
            type = unique(type),
            nContigs = dplyr::n(),
            totalLength = sum(length),
            nCircular = sum(circular == "Yes"),
            goodAssembly = case_when(
              nCircular >= 1 & nContigs <= 5 ~ "Yes",
              nCircular >= 1 & (nContigs - nCircular) <= 2  & nContigs < 8 ~ "Yes",
              TRUE ~ "No")) 

goodIsos <-  contigstats %>%
  filter(goodAssembly == "Yes", 
         !is.na(isolate),
         # !grepl("long", type, ignore.case = T)
         ) %>%
  distinct(isolate) %>%
  `$`(isolate)

# combine both and plot ---------------------------------------------------

stats <- bind_rows(nanostats, illustats %>% 
                     filter(type == "moretrimmed") %>%
                     ungroup() %>% 
                     select(all_of(colnames(nanostats))))

nanoplt <- bind_rows(nanotmplen) %>%
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
                        nanodir,".pdf")),
       nanoplt,
       height = 5.3, width = wd)

# save out stats in csv ---------------------------------------------------

all.stats <- contigstats %>%
  left_join(
    nanostats %>%
      select(
        isolate,
        nano.n = n,
        nano.N50 = N50,
        nano.mean = mean,
        nano.median = median,
        nano.total_bp = total_bp
      )
  ) %>%
  left_join(
    illustats %>%
      ungroup() %>%
      select(
        isolate,
        type,
        illu.n.contig = n,
        illu.N50 = N50,
        illu.mean = mean,
        illu.median = median,
        illu.total_bp = total_bp
      )
  ) %>%
  select(-guppy)

write_csv(all.stats, file.path(outdir,paste0("summary.", nanodir,".csv")))
