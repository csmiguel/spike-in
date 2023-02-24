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
# turn into proportions.
ps16s <-
  readRDS("data/intermediate/ps16s_ktu_filt_rel.rds") %>%
  phyloseq::transform_sample_counts(function(x) x / sum(x))
psITS <-
  readRDS("data/intermediate/psITS_ktu_filt_rel.rds") %>%
  phyloseq::transform_sample_counts(function(x) x / sum(x))
# remove spike-in sample
ab_ab <-
  abs_abundances %>%
  dplyr::filter(sample != "mcs168")

# NO CORRECTION for taxon-specific copy number is applied: it is assumed that the
# copy number of the marker is the same in the spike-in and species in the sample.
### 4. multiply rel otu table by total biological units in sample #####
# 16S
# multiply the counts in sample n by size_factors[n] and save in matrix `mat`
assertthat::assert_that(all(rownames(otu_table(ps16s)) == ab_ab$sample),
                        msg = "check names")
mat16s <-
  sweep(otu_table(ps16s), 1 + taxa_are_rows(ps16s), ab_ab$wild_16s_gsoil_av, FUN = "*")

# Make a new phyloseq object with the transformed counts
ps16s_abs <- ps16s
otu_table(ps16s_abs) <- otu_table(mat16s, taxa_are_rows = taxa_are_rows(ps16s)) %>% round()

# ITS
# multiply the counts in sample n by size_factors[n] and save in matrix `mat`
assertthat::assert_that(all(rownames(otu_table(psITS)) == ab_ab$sample),
                        msg = "check names")
matits <-
  sweep(otu_table(psITS), 1 + taxa_are_rows(psITS), ab_ab$wild_ITScell_gsoil, FUN = "*")
# Note, the result is a plain matrix and not an `otu_table` object

# Make a new phyloseq object with the transformed counts
psITS_abs <- psITS
otu_table(psITS_abs) <- otu_table(matits, taxa_are_rows = taxa_are_rows(psITS)) %>% round()

# save phyloseq object with corrected absolute counts
saveRDS(ps16s_abs, "data/intermediate/ps16s_abs.rds")
saveRDS(psITS_abs, "data/intermediate/psITS_abs.rds")
