# filter ktu phyloseq
library(phyloseq)
library(tidyverse)

#load phyloseq objects
load("data/intermediate/psKTU.Rdata")

#load functions for filtering phyloseq objects
dir("code/functions", "ps_filter", full.names = T) %>% sapply(source)
source("code/functions/otu_table2df.r")
source("code/functions/ps_merge.r")

# load ktu's corresponding to internal controls
ktu_spike_16s <- read.csv("output/spikein_summary_16s.csv")$OTU

# identify internal controls
int_controls <- "mcs168"

# 16S
sink(file = "output/filter_phyloseq_16s.txt")

# start filtering
ps16sktu_filt <-
  ps_16sktu %>%
  ps_filter_phylum_is_NA() %>%
  ps_filter_organelles() %>%
  # move internal controls to another ps, so they do not interfere with prevalence/abundance filtering
  ps_filter_split_taxa_ps(asv2second_ps = ktu_spike_16s, nameSecond_ps = "ps_16sktu_spike") %>%
  ps_filter_contaminants(blanks = "Bpcr|blank", excl = int_controls) %>%
  ps_filter_prevalence(
    mult_threshold = log(nsamples(.)) * 50, #I use log to account for
    # the fact that with larger sample size the chances of having shared false
    # positives across samples could increase in a logarithmic manner.
    prevalence_threshold = 2
  ) %>%
  ps_filter_relative_abundance(mean_prop = 5e-5) %>%
  # merge internal controls to filtered phyloseq
  {phyloseq(
    otu_table(merge_otutables(ps1 = ., ps2 = ps_16sktu_spike), taxa_are_rows = F),
    tax_table = merge_taxtable(ps_16sktu_spike, ., ps_16sktu),
    sample_data(ps_16sktu),
    refseq = c(refseq(ps_16sktu_spike), refseq(.))
  )}

sink()

# retain only samples [No internal controls] with no internal controls
ps16s_ktu_filt_rel <-
  subset_samples(ps16sktu_filt,
                !sample_names(ps16sktu_filt) %in% int_controls) %>%
                subset_taxa(!taxa_names(ps16sktu_filt) %in% ktu_spike_16s) %>%
                ps_filter_prevalence(prevalence_threshold = 2)

#save filtered phyloseq
saveRDS(ps16sktu_filt, "data/intermediate/ps16s_ktu_filt.rds")
saveRDS(ps16s_ktu_filt_rel, "data/intermediate/ps16s_ktu_filt_rel.rds")
