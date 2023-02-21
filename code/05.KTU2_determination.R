# KTU determination
# I will use beta KTU2 to cluster ASVs
library(tidyverse)
library(KTU2) #https://github.com/csmiguel/KTU2-beta-dada on commit 4b4d935
load("data/intermediate/seqtabNoC.Rdata")

# fungi
# since the sequence number is low I will run klustering instead of ktusp
fung <-
  dada2KTU(seqtab = seqtabNoC_ITS,
           cores = 2,
           method = "klustering"
  )
# bacteria
bact <-
  dada2KTU(seqtab = seqtabNoC_16S,
           cores = 2,
           split_lwrlim = 2000,
           split_reassemble = 500,
           method = "ktusp")

saveRDS(fung, "data/intermediate/ktu_fung.rds")

saveRDS(bact, "data/intermediate/ktu_bact.rds")