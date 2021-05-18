library(tidyverse)
library(yingtools2)
library(readxl)
library(gridExtra)
library(grid)

# argopts <- commandArgs(trailingOnly = TRUE)
# path = argopts[1]

dtstr <- gsub("^20","",str_extract(basename(getwd()), "^[0-9]+" ))

trimraw <- read_delim("trim.stats.txt",delim = ": ", 
                      col_names = c("logfn","time1","time2","time3",
                                    "readfoo","handlerfoo",
                                    "description","count"), trim_ws = T)

trimclean <- trimraw %>%
  transmute(logfn, countType = gsub("pair[12]","pair",handlerfoo), 
            description, count) %>%
  distinct(countType, count, .keep_all = T) %>%
  mutate(seq_id = sapply(str_split(logfn,"_"), function(x) x[1])) %>%
  # count(seq_id) %>% 
  select(seq_id, countType, count, description) 

cat("Saving out trimStats.rds for future postgres loading...\n")
saveRDS(trimclean,"trimStats.rds")

trimclean %>% 
  filter(grepl("[finlraw] pair", countType)) %>%
  select(-description) %>%
  spread(countType, count) %>%
  summarise_at(vars("final pair","raw pair"), 
               list(mean=mean, median = median))

plt <- trimclean %>%
  filter(grepl("pair",countType)) %>%
  mutate(countType = factor(countType, levels = c("raw pair",
                                                  "trimmed pair",
                                                  "decontaminated mouse_C57BL_6NJ pair",
                                                  "decontaminated hg37dec_v0.1 pair",
                                                  "final pair"),
                            labels = c("raw","trimmed","deMouse","deHuman","final"))) %>%
  filter(countType %in% c("raw","trimmed","deHuman","final")) %>%
  select(seq_id, countType, count) %>%
  spread(countType, count) %>%
  separate(seq_id, c("foodir","ID","fdir"), sep = "/") %>%
  mutate(ID = gsub("^[0-9]{6}-","", ID)) 

pgg <- plt %>%
  mutate(Trimmed = (raw - trimmed)/raw,
         Host = (trimmed - final)/raw,
         Clean = final/raw) %>%
  select(ID, Trimmed, Host, Clean) %>%
  gather( "type","perc",-ID) %>%
  ggplot(aes(ID, perc)) +
  geom_col(aes(fill = type)) +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  labs(x = "", y = "read %")

cgg <- plt %>%
  select(ID, raw, final) %>%
  gather("readType","cnt",-ID) %>%
  mutate(readType = factor(readType, levels = c("raw","final"))) %>%
  ggplot(aes(ID, cnt)) +
  geom_col(color = "black", fill = "grey", alpha = 0.65) +
  facet_grid(readType ~ .) +
  theme_bw() +
  # scale_y_continuous(trans = "log10") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Samples", y = "# of reads") +
  geom_hline(aes(yintercept = 1e6), color = "red", linetype = 2) +
  geom_hline(aes(yintercept = 5e6), color = "blue", linetype = 2)

pdf(paste0(dtstr,".readqual.pdf"), height = 7.2, width = 11)

gg.stack(pgg, cgg, heights = c(4,3),newpage = F)

dev.off()
