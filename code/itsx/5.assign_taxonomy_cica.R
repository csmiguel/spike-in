# assign taxonomy
library(tidyverse)
library(Biostrings)
library(DECIPHER)

# ktu object (list) to input for idTaxa (matrix)
source("/home/mcamacho/spike-in/klustering2otutable.r")
# load customized idTaxa
source("/home/mcamacho/spike-in/assign_taxonomy_decipher.R")
# load function to clean ITS output
source("/home/mcamacho/spike-in/clean_idTaxa_its.R")

# load KTU its2
ktuits <-
  readRDS("/home/mcamacho/spike-in/ktu_its2.rds") %>%
  klustering2otutable()

# assign taxonomy
#load UNITE trainingSet
load("/home/mcamacho/spike-in/UNITE_v2021_May2021.RData")
# assign taxonomy
taxid_its <-
  assign_taxonomy_decipher(seqtabdada2 = ktuits,
                           training_set = trainingSet,
                           nprocessors = 1) %>%
  clean_its_idTaxa()
rm(trainingSet)
#taxid objects have the below structure. rownames are seqs and taxonomy in cols.
# (rownames) domain  phylum    class      order     family    genus     species
# ACCTAT… Bacter… Campilobac… Campylob… Campylo… Sulfurovace… Sulfurovum    NA
# ACCTCT… Bacter… Proteobact… Alphapro… Defluvi… Defluviicoc… Defluviicocc… NA
# ACCTGT… Bacter… Desulfobac… Desulfob… Desulfo… Desulfosarc… NA            NA

#save seqs with assigned taxonomy
saveRDS(taxid_its, "/home/mcamacho/spike-in/taxid_its2.rds")
