# generate nanopore long and Illumina reads length distribution and hybrid assembly results

library(stringr)
library(tidyverse)
library(lubridate)
library(conflicted)

argoptions <- commandArgs(trailingOnly = TRUE)
nanopath <- argoptions[1]
nanopath <- gsub("/$","",nanopath)
nanodir <- basename(nanopath)

if (is.na(argoptions[3])){
  outdir <- getwd()
} else {
  outdir <- argoptions[3]
}

if (is.na(argoptions[2])){
  name.idx <- 9
} else {
  name.idx <- 10
}

conflict_prefer("n", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")

calNx <- function(lenArr, nx = 0.5){
  len.sorted <- rev(sort(lenArr))
  len.sorted[cumsum(len.sorted) >= sum(len.sorted)*nx][1]
}

calWidth <- function(nsamples){ return(5 + 1.3*nsamples) }

# nanopath <- "/Volumes/dfi-cores/DFI-MMF/Nanopore/nanopore/Bioinformatics/20240207_ONT67"

# stats for nanopore long reads ------------------------------------------------
nanofais <- Sys.glob(file.path(nanopath,"*/*.fastq.fai"))
# nanofais <- nanofais[grepl("barcode[0-9]{2}.fastq.fai",nanofais)]

nanotmplen <- NULL

cat(">Reading Nanopore reads indexes....\n")

for (f in nanofais){

  iso <- str_split(f,"/")[[1]][name.idx]

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

# nanostats in tsv for long reads -----------------------------------------

nstatfais <- Sys.glob(file.path(nanopath,"*/NanoStat.output"))

cat(">Reading NanoStat output....\n")

nstat <- lapply(nstatfais, function(fn) {read_tsv(fn) %>% mutate(filepath = fn,
                                                        isolate = str_split(fn, "/")[[1]][name.idx])}) %>% 
  bind_rows()

nstat.org <- nstat %>% 
  filter(Metrics %in% c("read_length_stdev",
                        "mean_qual",
                        "mean_read_length",
                        "median_qual",
                        "median_read_length",
                        "number_of_reads")) %>% 
  transmute(Metrics, isolate, dataset = as.numeric(dataset)) %>% 
  spread(Metrics, dataset)

# Illumina -----------------------------------------------------

illufais <- Sys.glob(file.path(nanopath,"*/*/002*.fai"))
# illufais <- illufais[grepl("hybridAssemblyTrim",illufais)]

illutmplen <- NULL

cat(">Reading Illumina contigs indexes....\n")

for (f in illufais){

  iso <- str_split(f,"/")[[1]][name.idx]

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

cat(">Reading flye results....\n")

for (flyfn in flye.fns){
  
  iso <- str_split(flyfn,"/")[[1]][name.idx]
  
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

# get pilon location if there is any --------------------------------------

pilonfais <- Sys.glob(file.path(nanopath,"*/*/pilon.fasta.fai"))

pilontmplen <- NULL

cat(">Reading pilon results....if any....\n")

for (f in pilonfais){
  
  iso <- str_split(f,"/")[[1]][name.idx]
  
  pilontmplen[[iso]] <- read_tsv(f, col_names = F) %>%
    transmute(contigNo= gsub("_pilon","",X1), 
              pilon.length = X2,
              pilon.path = gsub(".fai$","",f)) %>%
    mutate(isolate = iso)
  
}

pilon.len <- bind_rows(pilontmplen) %>% 
  mutate(pilon.path = gsub("/Volumes/dfi-cores/DFI-MMF/Nanopore/nanopore/Bioinformatics","## .", pilon.path))

# summary of the assemblies -----------------------------------------
sumfn <- Sys.glob(file.path(nanopath,"summary*txt"))

cat(paste0(">Reading in the latest summary file: ",last(sumfn),".....\n"))
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
  relocate(isolate, type) %>% 
  left_join(pilon.len)

cat(paste0(">>Output detailed summary file: ",
           paste0("detailed.summary.", nanodir,".nano.csv"),
           ".....\n"))

write_csv(org.summ %>% 
            select(-guppy), 
          file.path(outdir,paste0("detailed.summary.", nanodir,".nano.csv")))

contigstats <- org.summ %>%
  mutate(int.length = if_else(is.na(pilon.path),
                              length,
                              pilon.length),
         int.path = if_else(is.na(pilon.path),
                            path,
                            pilon.path)) %>% # view
  group_by(int.path) %>%
  summarise(guppy = unique(guppy),
            isolate = unique(isolate),
            type = unique(type),
            nContigs = dplyr::n(),
            totalLength = unique(sum(int.length)),
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

cat(paste0(">>Output plot: ",
           paste0("nanopore.isolates.hist.",
                  nanodir,".pdf"),
           ".....\n"))

ggsave(file.path(outdir,
                 paste0("nanopore.isolates.hist.",
                        nanodir,".pdf")),
       nanoplt,
       height = 5.3, width = wd)

# save out stats in csv ---------------------------------------------------

all.stats <- contigstats %>%
  dplyr::rename(path = int.path) %>% 
  left_join(
    nstat.org
  ) %>% 
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
  select(-guppy) %>% 
  arrange(isolate)

cat(paste0(">>Output final summary file: ",
           paste0("summary.", nanodir,".csv"),
           ".....\n"))

write_csv(all.stats, file.path(outdir,paste0("summary.", nanodir,".csv")))
