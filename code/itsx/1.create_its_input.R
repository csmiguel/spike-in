library(tidyverse)
library(Biostrings)
library(digest)

load("data/intermediate/seqtabNoC.Rdata")
itsx <- "data/intermediate/ITSx"

if(!dir.exists(itsx))
   dir.create(itsx)

# ktu fung is the product of applying KTU2 to the seqtable
its_seqs <-
  seqtabNoC_ITS %>%
  colnames() %>%
  Biostrings::DNAStringSet()
names(its_seqs) <-
  its_seqs %>%
      as.character %>%
      sapply(digest, "md5")
# write
Biostrings::writeXStringSet(
  x = its_seqs,
  filepath = file.path(itsx, "input_itsx.fasta"),
  format = "fasta"
)
