# produce seqtab from ITSx output so as to stream it into the general workflow
library(tidyverse)
library(Biostrings)
library(digest)

# read extracted ITS2
its2 <- 
  Biostrings::readDNAStringSet("data/intermediate/ITSx/its2_itx.fasta")
names(its2) <-
  names(its2) %>% stringr::str_remove("\\|.*$")

load("data/intermediate/seqtabNoC.Rdata")

# if not digested
# digest seqtab DNAseq to md5.
if(nchar(colnames(seqtabNoC_ITS)[1]) > 100) {
  colnames(seqtabNoC_ITS) <-
      colnames(seqtabNoC_ITS) %>%
      sapply(digest, "md5") %>% {
        attributes(.) <- NULL; .
      }
  }

# filter ITSx fungi sequences from the seqtab object
# keep only those present in the extracted ITS2 sequences.
# replace the names in the matrix to the ITS2 extracted sequence.
seqtab_its2 <-
  seqtabNoC_ITS[, names(its2)] %>%
  {
    colnames(.) <-
      its2[colnames(.)] %>%
      as.character %>% {
      attributes(.) <- NULL; .
      } ; .
  }

# save
saveRDS(seqtab_its2, "data/intermediate/ITSx/seqtab_its2.rds")
