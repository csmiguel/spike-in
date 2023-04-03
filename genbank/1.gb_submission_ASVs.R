library(phyloseq)
library(tidyverse)
library(ShortRead)
library(xlsx)

### 16S ####
### 1. write FASTA files ###
# bact
ps_bact <-
  readRDS("data/intermediate/ps16s_abs.rds")

refseq(ps_bact) %>%
  ShortRead::writeFasta("genbank/16s.fasta")
# read biosample info
bios <-
  xlsx::read.xlsx("data/raw/biosamples.xlsx", 1) %>%
  dplyr::select(biosample_accession, sample_name) %>%
  dplyr::filter(sample_name %in% sample_names(ps_bact))
### 2. create mapping files ###
# 16s
# chimeric ASVs detected by NCBI submission wizard
seqs2remove <-
  readLines("genbank/SUB13022362_chimera_scores.txt") %>%
  {.[grep("ASV", .)]} %>%
  stringr::str_replace("\t.*", "")
  
# format metadata table for genbank
mapping_complete <-
  ps_bact %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  dplyr::select(OTU, Sample) %>%
  left_join(bios, by = c("Sample" = "sample_name")) %>%
  dplyr::select(-Sample) %>%
  plyr::ddply(~OTU, function(x) {
    data.frame(Sequence_ID = unique(x$OTU),
               biosample_accession = paste(x$biosample_accession, collapse = ","))
  }) %>%
  select(Sequence_ID, biosample_accession)  %>%
  dplyr::filter(!Sequence_ID %in% seqs2remove)

# the submission wizard only accepts one biosample per ASV
mapping_1sample <-
  mapping_complete %>%
  mutate(biosample_accession = stringr::str_remove(biosample_accession, ",.*$"))
# write files
write_tsv(mapping_complete, "genbank/mapping_file16s_complete.tsv")
write_tsv(mapping_1sample, "genbank/mapping_file16s_oneSample.tsv")


### ITS ####

### 1. write FASTA files ###
# fungi  
ps_fung <-
  readRDS("data/intermediate/psITS_abs.rds")

refseq(ps_fung) %>%
  ShortRead::writeFasta("genbank/its.fasta")

### 2. create mapping files ###
# format metadata table for genbank
mapping_complete_its <-
  ps_fung %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  dplyr::select(OTU, Sample) %>%
  left_join(bios, by = c("Sample" = "sample_name")) %>%
  dplyr::select(-Sample) %>%
  plyr::ddply(~OTU, function(x) {
    data.frame(Sequence_ID = unique(x$OTU),
               biosample_accession = paste(x$biosample_accession, collapse = ","))
  }) %>%
  select(Sequence_ID, biosample_accession)

# write files
write_tsv(mapping_complete_its, "genbank/mapping_fileITS_complete.tsv")
