# KTU determination on ITS2 sequences
# I will use beta KTU2 to cluster ASVs
library(tidyverse)
library(KTU2) #https://github.com/csmiguel/KTU2-beta-dada on commit 4b4d935
seqtab <- readRDS("data/intermediate/ITSx/seqtab_its2.rds")

# since the sequence number is low I will run klustering instead of ktusp
fung <-
  dada2KTU(seqtab = seqtab,
           cores = 2,
           method = "klustering"
  )
# save ktu object
saveRDS(fung, "data/intermediate/ITSx/ktu_its2.rds")
