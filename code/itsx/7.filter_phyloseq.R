# detect spike-in
library(tidyverse)
library(phyloseq)

ps_its2 <- readRDS("data/intermediate/ITSx/psKTU.rds")
source("code/functions/ps_filter_prevalence.r")

# spike-in ----------------------------------------------------------------
# detect spike-in
spikeits <-
  ps_its2 %>%
  prune_samples(samples = "mcs168") %>%
  ps_filter_prevalence(abundance_threshold = 1) %>%
  tax_glom("family") %>%
  # yarrowia was not detected
  subset_taxa(family %in% "Dipodascaceae") %>%
  psmelt()


# save results
write.csv(spikeits,
          file = "output/itsx_spikein_summary_its.csv",
          row.names = F)

# filter phyloseq ----------------------------------------------------------------
#load functions for filtering phyloseq objects
dir("code/functions", "ps_filter", full.names = T) %>% sapply(source)
source("code/functions/otu_table2df.r")
source("code/functions/ps_merge.r")

# load ktu's corresponding to internal controls
ktu_spike_its <- read.csv("output/itsx_spikein_summary_its.csv")$OTU

# identify internal controls
int_controls <- "mcs168"

#open connection
sink(file = "output/ITSx_filter_phyloseq_ITS.txt")
# start filtering
psITS_ktu_filt <-
  ps_its2 %>%
#  ps_filter_phylum_is_NA() %>%
  # move internal controls to another ps, so they do not interfere with prevalence/abundance filtering
  ps_filter_split_taxa_ps(asv2second_ps = ktu_spike_its, nameSecond_ps = "ps_itsktu_spike") %>%
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
    otu_table(merge_otutables(ps1 = ., ps2 = ps_itsktu_spike), taxa_are_rows = F),
    tax_table = merge_taxtable(ps_itsktu_spike, ., ps_its2),
    sample_data(ps_its2),
    refseq = c(refseq(ps_itsktu_spike), refseq(.))
  )}

sink()

# retain only samples [No internal controls]
psITS_ktu_filt_rel <-
  subset_samples(psITS_ktu_filt,
                 !sample_names(psITS_ktu_filt) %in% int_controls) %>%
  subset_taxa(!taxa_names(psITS_ktu_filt) %in% ktu_spike_its) %>%
  ps_filter_prevalence(prevalence_threshold = 2)

#save filtered phyloseq
saveRDS(psITS_ktu_filt, "data/intermediate/ITSx/psITS2_ktu_filt.rds")
saveRDS(psITS_ktu_filt_rel, "data/intermediate/ITSx/psITS2_ktu_filt_rel.rds")
