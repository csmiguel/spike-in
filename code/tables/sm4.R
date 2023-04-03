# add sample info to genbank accession table
library(tidyverse)
library(xlsx)

# biosample data
bios <-
  xlsx::read.xlsx("data/raw/biosamples.xlsx", 1)
# sample data
samples <-
  readRDS("data/intermediate/samples.rds")
# metadata ITS SRA
sraits <-
  read.csv("data/raw/metadata_SRA_its.tsv", sep = "\t") %>%
  rename(library = library_ID,
         SRA_ITS2_accession = accession) %>%
  select(library, SRA_ITS2_accession)

# merge
bios %>%
  left_join(samples,
            by = "code") %>%
  dplyr::rename(library = sample_name,
                SRA_16S_accession = accession) %>%
  dplyr::left_join(sraits,
                   by = "library") %>%
  dplyr::select(sample, library, code,
                bioproject_accession, biosample_accession,
                SRA_16S_accession, SRA_ITS2_accession) %>%
    write.xlsx("output/sm4.xlsx")
