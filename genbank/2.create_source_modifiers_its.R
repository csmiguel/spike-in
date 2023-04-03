library(phyloseq)
library(tidyverse)
library(ShortRead)

#phyloseq taxonomy
# fungi  
ps_fung <-
  readRDS("data/intermediate/psITS_abs.rds")

taxt <-
  tax_table(ps_fung) %>%
  as.data.frame %>%
  select(-species)

taxt[apply(taxt, 2, function(x) grepl(pattern = "unidentified", x = x))] <- NA

# get organism names
organisms <-
  taxt %>%
  apply(1, function(x) {
    last_rank <-
      x[length(x) - sum(is.na(x))] %>%
      sub(pattern = "_.*$", replacement = "")
    paste("uncultured", last_rank)
  })

# submitted fasta
fungreads <-
  ShortRead::readFasta("genbank/its.fasta")
# assert all fasta seqs are in ps
if(!all(sort(as.character(ShortRead::id(fungreads))) == sort(names(organisms))))
  error("fasta seqs are not the same than those in the phyloseq object")

# build data frame for source modifier
data.frame(
  Sequence_ID = names(organisms),
  Organism = organisms,
  clone = names(organisms),
  isolation_source = "soil from strawberry plantations"
  ) %>%
  write_tsv("genbank/source_modifiers_its.tsv")
