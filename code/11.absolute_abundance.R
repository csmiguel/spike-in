# compute absolute abundances
# output is a phyloseq object with absolute sample counts
library(phyloseq)
library(plyr)
library(tidyverse)

source("code/functions/absolute_ab.R")

# load filtered phyloseq with no blanks but all the internal controls
ps16s <-
  readRDS("data/intermediate/ps16s_ktu_filt.rds")
psITS <-
  readRDS("data/intermediate/psITS_ktu_filt.rds")

# ASV internal controls
ktu_spike_16s <- read.csv("output/spikein_summary_16s.csv")
ktu_spike_its <- read.csv("output/spikein_summary_its.csv")

###1. estimate spike-in marker copies in inoculated sample#####
# stock concentrations
# copies of 16s per 20 ul of stock (units: cp/ul)
imte_c1 <- 6*10^7 / 20 # units: copies/ul
allo_c1 <- 1.4*10^8 / 20
# for the yeast Yarrowia I used the stock of 10^7 cells/ml
# 10^7 cells/ml * 1ml/1000ul
yarr_c1 <- 10^7 / 1000 # units = cells/ul

# concentration in spike-in (c2)
# volume of zymo research stock in spike in solution
v1_zymo <- 265
# volume of yarrowia stock 10^7 in spike-in solution
v1_yarrowia <- 662
# total volume of spike in solution
v2_spike <- v1_zymo + v1_yarrowia

# concentration of spiked-in bact/yeast in spike in solution [16s copies]|[cells] per ul.
allo_c2 <- c2(c1 = allo_c1, v1 = v1_zymo, v2 = v2_spike)
imte_c2 <- c2(c1 = imte_c1, v1 = v1_zymo, v2 = v2_spike)
yarr_c2 <- c2(c1 = yarr_c1, v1 = v1_yarrowia, v2 = v2_spike)

# Kc, spike-in marker copies in inoculated sample.
Kc <-
  sample_data(ps16s) %>%
  data.frame() %>%
  mutate(
    allo_16s = spikein * allo_c2,
    imte_16s = spikein * imte_c2,
    yarrowia_cel = spikein * yarr_c2) %>%
  tibble::rownames_to_column("sample")

###2. estimate wild copies in samples#####
abs_abundances <-
  # reads per spike in and total
  list(total_reads_16s =
         otu_table(ps16s) %>% rowSums(), # total reads per sample
       total_reads_its =
         otu_table(psITS) %>% rowSums(), # total reads per sample
       allo_reads =
         prune_taxa(ktu_spike_16s$OTU[ktu_spike_16s$genus == "Allobacillus"], ps16s) %>%
         otu_table() %>% rowSums(),# total reads per sample of Allobacillus
       imte_reads =
         prune_taxa(ktu_spike_16s$OTU[ktu_spike_16s$genus == "Imtechella"], ps16s) %>%
         otu_table() %>% rowSums(), # total reads per sample of Imtechella
       yarrowia_reads =
         prune_taxa(ktu_spike_its$OTU[1], psITS) %>%
         otu_table() %>% rowSums() # total reads per sample Yarrowia
  ) %>%
  # I am obliged to do this trick to merge the results from the list into a DF.
  lapply(function(x) data.frame(sample = names(x), reads = x)) %>%
  plyr::ldply(.id = "var") %>%
  pivot_wider(names_from = var, values_from = reads) %>%
  # lastly, I compute the total number of wild reads. (wild == non spike in)
  mutate(wild_reads_16s = total_reads_16s - allo_reads - imte_reads,
         wild_reads_its = total_reads_its - yarrowia_reads) %>%
  left_join(Kc, by = "sample") %>%
  # compute copies off wild (non-control) 16s copies/ITS equiv. cells in 1 gram of soil.
  mutate(
    wild_16s_gsoilAllo = wild_copies_gsoil(sample_wild_reads = wild_reads_16s, # number of wild 16S copies in the soil following calculations with Allobacillus
                                           spikein_reads = allo_reads,
                                           spikein_copies = allo_16s,
                                           soilg = soil_mg / 1000),
    wild_16s_gsoilImt = wild_copies_gsoil(sample_wild_reads = wild_reads_16s, # number of wild 16S copies in the soil following calculations with Imtechella
                                          spikein_reads = imte_reads,
                                          spikein_copies = imte_16s,
                                          soilg = soil_mg / 1000),
    wild_ITScell_gsoil = wild_copies_gsoil(sample_wild_reads = wild_reads_its, # number of yarrowia-equivalents ITS copies in the soil.
                                           spikein_reads = yarrowia_reads,
                                           spikein_copies = yarrowia_cel,
                                           soilg = soil_mg / 1000),
    wild_16s_gsoil_av = (wild_16s_gsoilAllo + wild_16s_gsoilImt) / 2 # average wild 16s copies from Imtechella and Allobacillus
  ) %>%
  as_tibble()

saveRDS(abs_abundances, "data/intermediate/abs_abundances.rds")

###3. estimate scaling factors #####

rm(ps16s); rm(psITS)
# load filtered phyloseq with no blanks and no internal controls
ps16s <-
  readRDS("data/intermediate/ps16s_ktu_filt_rel.rds")
psITS <-
  readRDS("data/intermediate/psITS_ktu_filt_rel.rds")
# remove spike-in sample
ab_ab <-
  abs_abundances %>%
  dplyr::filter(sample != "mcs168")

# NO CORRECTION for taxon-specific copy number is applied
# 16S
# create vector with size factors and order as sample names
# the result after correction will be the number of 16 copies in 1g of soil
size_factors16s_gsoil <- {
  m <- match(sample_names(ps16s), ab_ab$sample)
  # ratio copies in soil sample respect to reads
  r <- ab_ab$wild_16s_gsoil_av[m] / ab_ab$wild_reads_16s[m]
}

# ITS
# create vector with size factors and order as sample names
# wild_copies / wild_reads: the result after correction will be the number of ITS-yarrowia equivalent copies in 1g of soil.
size_factorsITS_gsoil <- {
  m <- match(sample_names(psITS), ab_ab$sample)
  # ratio copies in soil sample respect to reads
  r <- ab_ab$wild_ITScell_gsoil[m] / ab_ab$wild_reads_its[m]
}

### 4. muliply otu table by scaling factors #####
# 16S
# multiply the counts in sample n by size_factors[n] and save in matrix `mat`
mat16s <-
  sweep(otu_table(ps16s), 1 + taxa_are_rows(ps16s), size_factors16s_gsoil, FUN = "*")
# Note, the result is a plain matrix and not an `otu_table` object

# Make a new phyloseq object with the transformed counts
ps16s_abs <- ps16s
otu_table(ps16s_abs) <- otu_table(mat16s, taxa_are_rows = taxa_are_rows(ps16s))

# ITS
# multiply the counts in sample n by size_factors[n] and save in matrix `mat`
matits <-
  sweep(otu_table(psITS), 1 + taxa_are_rows(psITS), size_factorsITS_gsoil, FUN = "*")
# Note, the result is a plain matrix and not an `otu_table` object

# Make a new phyloseq object with the transformed counts
psITS_abs <- psITS
otu_table(psITS_abs) <- otu_table(matits, taxa_are_rows = taxa_are_rows(psITS))

# save phyloseq object with corrected absolute counts
saveRDS(ps16s_abs, "data/intermediate/ps16s_abs.rds")
saveRDS(psITS_abs, "data/intermediate/psITS_abs.rds")

