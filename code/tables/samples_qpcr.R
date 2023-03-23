# table with qpcr samples
library(tidyverse)
idqpcr <-
  readRDS("data/intermediate/id_qpcr6.rds") %>%
  arrange(code) %>%
  rename(dna_ext = sample)
sampleid <- readRDS("data/intermediate/samples.rds")
sampleid %>%
  left_join(idqpcr,
            by = "code") %>%
  drop_na() %>%
  write.table("output/samples_qpcr.txt")

