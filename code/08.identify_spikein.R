library(tidyverse)
library(phyloseq)

load("data/intermediate/psKTU.Rdata")
source("code/functions/ps_filter_prevalence.r")

# read and filter phyloseq
# I use a trick to find internal controls from the mock community, consisting in keeping only both Zymo samples
# and retaining those KTU present in both samples and with a least 1 read.
# 16S
spike16s <-
  ps_16sktu %>%
  prune_samples(samples = "mcs168") %>%
  ps_filter_prevalence(abundance_threshold = 1) %>%
  tax_glom("genus") %>%
  subset_taxa(genus %in% c("Imtechella", "Allobacillus")) %>%
  psmelt()

# ITS
spikeits <-
  ps_itsktu %>%
  prune_samples(samples = "mcs168") %>%
  ps_filter_prevalence(abundance_threshold = 1) %>%
  tax_glom("genus") %>%
  psmelt()

# save results
write.csv(spikeits,
          file = "output/spikein_summary_its.csv",
          row.names = F)
write.csv(spike16s,
          file = "output/spikein_summary_16s.csv",
          row.names = F)
