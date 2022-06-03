setwd("~/OneDrive - The University of Chicago/16s/libQual")

# 2020.7.7
# check high quality survival percentage across ALL libraries

rm(list = ls())

library(yingtools2)
library(readr)
library(tidyverse)
library(phyloseq)

# addIDs <- function(file){
#   tmpCSV <- read_csv(file)
#   pathSplit <- strsplit(file, "/")[[1]]
#   pool <- pathSplit[length(pathSplit)-2]
#   tmpCSV <- tmpCSV %>%
#     mutate(sampleid = paste0(pool, "_", gsub("_R1.fastq","", X1)),
#            pool = pool)
#   return(tmpCSV)
# }

# GET ALL FILTER STATS
# filt_files <- Sys.glob(file.path("/Volumes/pamer-lab/DFI_MMF/MiSeq/Bioinformatics/runPools",
#                                    "*/SM*/cleanFastq/filterTrim.stats.csv"))
# outfilt_files <- Sys.glob(file.path("/Volumes/dfi-cores/DFI-MMF/MiSeq/projects/*",
#                    "SM*/cleanFastq/filterTrim.stats.csv"))
#
# fil_list <- lapply(c(filt_files,outfilt_files), addIDs)

fiphy <- readRDS("/Volumes/pamer-lab/DFI_MMF/MiSeq/Bioinformatics/dada2rds/finalPhy_rdp.merged.wMeta.rds")

# IF NEED TO FIX META FILES
# meta <- readRDS("/Volumes/pamer-lab/DFI_MMF/MiSeq/Bioinformatics/dada2rds/meta.merged.rds")
#
# cmeta <- meta %>%
#   as.matrix() %>%
#   as.data.frame() %>%
#   mutate_if(is.factor,as.character)
#
# cmeta <- type_convert(as_tibble(cmeta), col_types = COLTYPES)
# newmeta <- sample_data(cmeta)
# sample_names(newmeta) <- newmeta$samplename
#
# newmeta <- sample_data(meta) %>%
#   as.matrix() %>%
#   as.data.frame() %>%
#   dplyr::select(-reads.in, -reads.out) %>%
#   left_join(bind_rows(fil_list) %>%
#               dplyr::rename(libraryname = pool,
#                             samplename = sampleid)) %>%
#   dplyr::select(-X1)
#
# rownames(newmeta) <- newmeta$samplename
# newmeta <- sample_data(newmeta)
#
# sample_data(fiphy) <- newmeta
# saveRDS(newmeta, "/Volumes/pamer-lab/DFI_MMF/MiSeq/Bioinformatics/dada2rds/meta.merged.rds")
# # sample_names(newmeta)
# saveRDS(fiphy, "/Volumes/pamer-lab/DFI_MMF/MiSeq/Bioinformatics/dada2rds/finalPhy_rdp.merged.wMeta.rds")
# seqtab <- readRDS("/Volumes/pamer-lab/DFI_MMF/MiSeq/Bioinformatics/dada2rds/seqtab.merged.rds")

totals <- sample_sums(fiphy) %>%
  as.data.frame() %>%
  rownames_to_column(var = "samplename") %>%
  as_tibble() %>%
  dplyr::rename(after = 2)

saN <- sample_data(fiphy) %>%
  data.frame() %>%
  as_tibble() %>%
  dplyr::select(samplename, libraryname) %>%
  dplyr::count(libraryname, name = "sampleN") %>%
  mutate(liborder = if_else(grepl("^SM",libraryname),
                            as.numeric(str_extract(libraryname,"[0-9]+")),
                            as.numeric(str_extract(libraryname,"[0-9]+$"))))

tmpStatsall <- sample_data(fiphy) %>%
  data.frame() %>%
  as_tibble() %>%
  dplyr::select(samplename, libraryname, reads.in, reads.out) %>%
  left_join(totals) %>%
  mutate(highQualSurvPerc = after/reads.in*100,
         filterSurPerc = reads.out/reads.in*100,
         liborder = if_else(grepl("^SM",libraryname),
                            as.numeric(str_extract(libraryname,"[0-9]+")),
                            as.numeric(str_extract(libraryname,"[0-9]+$"))))

pltdf <- tmpStatsall %>%
  dplyr::select(-after, -reads.out) %>%
  gather("precentage", "value", -samplename, -libraryname, -reads.in,-liborder)

pergg <- pltdf %>%
  ggplot(aes(x = reorder(libraryname,liborder), y = value)) +
  geom_boxplot(aes(fill = precentage),
               outlier.shape = 21, outlier.fill = "green",
               outlier.alpha = 0.55, outlier.size = 0.85,
               # position = position_dodge(0.95),
               ) +
  # stat_summary(aes(group = precentage),
  #              fun="mean", color="orange",
  #              geom = "line", linetype = 2,
  #              position = position_dodge(0.95))+
  theme_light() +
  labs(x = "Pool",
       y = "Survived %") +
  theme(axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_brewer(palette = "Set1")

totgg <- tmpStatsall %>%
  dplyr::select(samplename, libraryname, reads.in, after, liborder) %>%
  gather("precentage", "value", -samplename, -libraryname, -liborder) %>%
  mutate(precentage = factor(precentage, levels = c("reads.in","after"),
                             labels = c("raw","final high quality"))) %>%
  ggplot(aes(x = reorder(libraryname,liborder))) +
  geom_boxplot(aes(y = value, fill = precentage),
               outlier.shape = 21, outlier.fill = "green",
               outlier.alpha = 0.55, outlier.size = 0.85,) +
  geom_hline(yintercept = 5000, linetype = 2, color = "blue") +
  geom_hline(yintercept = 1000, linetype = 2, color = "red") +
  theme_bw() +
  labs(x = "Pool", fill = "Read type",
       y = "Total Raw/High quality\nReads/Sample") +
  scale_y_log10() +
  theme(axis.text.x.bottom = element_blank()) +
  paletteer::scale_fill_paletteer_d("suffrager::classic")

# totgg

ngg <- saN %>%
  ggplot(aes(x = reorder(libraryname,liborder))) +
  geom_col(aes(y = sampleN), color = "black",
           alpha = 0.65, width = 0.65, fill = "grey") +
  theme_bw() +
  theme(axis.text.x.bottom = element_blank()) +
  labs(x = " ",
       y = "# of Samples")

library(data.table)
library(gtable)
library(gridExtra)
library(grid)

wd <- 2.2 + 0.28* nrow(saN)

pdf(paste0("highqual_by_lib.",gsub("-","",lubridate::today()),".pdf"),
    height = 7.6,
    width = wd)
gg.stack(totgg, ngg, pergg, heights=c(2,1,3),newpage = F)
dev.off()

### COL TYPES HARD-CODED
COLTYPES <- cols(
  request.by = col_character(),
  lab.affiliation = col_character(),
  experiment.sample.name = col_character(),
  group = col_character(),
  mouse.number = col_character(),
  sample.weight = col_double(),
  treatment = col_character(),
  sample.location=col_character() ,
  species=col_character() ,
  mouse.strain=col_character() ,
  sample.type=col_character() ,
  sample.submission.date= col_datetime(),
  sampleid=col_character() ,
  forwardbarcode=col_character() ,
  reversebarcode=col_character() ,
  libraryname=col_character() ,
  adapter=col_character() ,
  pcr = col_character() ,
  miseq.submission= col_date(),
  miseq.data.receive.date= col_date(),
  reads.in = col_integer(),
  reads.out = col_integer(),
  samplename =col_character()
)
