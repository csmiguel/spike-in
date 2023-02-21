# assign taxonomy
library(tidyverse)
library(Biostrings)
library(DECIPHER)

# ktu object (list) to input for idTaxa (matrix)
source("code/functions/klustering2otutable.r")
# ktu objects
ktu16s <-
  readRDS("data/intermediate/ktu_bact.rds") %>%
  klustering2otutable()
  
ktuits <-
  readRDS("data/intermediate/ktu_fung.rds") %>%
  klustering2otutable()
# load customized idTaxa
source("code/functions/assign_taxonomy_decipher.R")
# load function to clean ITS output
source("code/functions/clean_idTaxa_its.R")

# assign taxonomy
load("data/raw/SILVA_SSU_r138_2019.RData")
taxid_16S <-
  assign_taxonomy_decipher(seqtabdada2 = ktu16s,
                           training_set = trainingSet,
                           nprocessors = 1)
rm(trainingSet) #I will remove it since both are named similarly

# fungi
#load UNITE trainingSet
load("data/raw/UNITE_v2020_February2020.RData")
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
save(taxid_16S, taxid_its, file = "data/intermediate/taxid.Rdata")
