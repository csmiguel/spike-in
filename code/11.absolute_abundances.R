# compute absolute abundances from phyloseq
# output is a sample_data data frame with extra columns with absolute abundances per sample
library(phyloseq)
library(tidyverse)
library(reshape2)

# load filtered phyloseq with no blanks but all the internal controls
ps16s <-
  readRDS("data/intermediate/ps16s_ktu_filt.rds")
psITS <-
  readRDS("data/intermediate/psITS_ktu_filt.rds")

# ASV internal controls
ktu_spike_16s <- read.csv("output/spikein_summary_16s.csv")
ktu_spike_its <- read.csv("output/spikein_summary_its.csv")

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

##functions###
# function to compute c2: c1.v1 = c2.v2
c2 <- function(c1, v1, v2) {
  h <- c1*v1/v2
  return(h)
}
# transform counts to reflect absolute abundances. Compute the total number of 16s copies/cells per gram of soil based on the spike.
wild_copies_gsoil <- function(sample_wild_reads, spikein_reads, spikein_copies, soilg) {
  # sample_wild_reads, non-spike in taxa in reads from the phyloseq object
  # spikein_reads, spike in reads in the phyloseq object
  # spikein_copies, number of spikein copies added to the soil sample at DNA extraction
  # soilg, grams of soil in DNA extraction
  sample_wild_copies = sample_wild_reads / spikein_reads * spikein_copies
  # for 1 gram of soil
  per_g_soil <- sample_wild_copies / soilg
  return(per_g_soil)
}
###

# concentration of spiked-in bact/yeast in spike in solution [16s copies]|[cells] per ul.
allo_c2 <- c2(c1 = allo_c1, v1 = v1_zymo, v2 = v2_spike)
imte_c2 <- c2(c1 = imte_c1, v1 = v1_zymo, v2 = v2_spike)
yarr_c2 <- c2(c1 = yarr_c1, v1 = v1_yarrowia, v2 = v2_spike)

# added units of spike-in to soil samples
qt <-
  sample_data(ps16s) %>%
  data.frame() %>%
  mutate(
    allo_16s = spikein * allo_c2,
    imte_16s = spikein * imte_c2,
    yarrowia_cel = spikein * yarr_c2) %>%
  tibble::rownames_to_column("sample")

# add read count and compute ITS/16s unit count on soil sample
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
  melt() %>%
  pivot_wider(names_from = L1, values_from = value) %>%
  # lastly, I compute the total number of wild reads. (wild == non spike in)
  mutate(wild_reads_16s = total_reads_16s - allo_reads - imte_reads,
         wild_reads_its = total_reads_its - yarrowia_reads) %>%
  left_join(qt, by = "sample") %>%
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
