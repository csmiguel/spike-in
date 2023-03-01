# plot expected versus observed proportion of spikein
library(tidyverse)
source("code/functions/kt_expected.R")

###1. observed proportion of spikein###
# proportion of observed spikein reads respect to total
obs <-
  readRDS("data/intermediate/abs_abundances.rds") %>%
  select(sample, code, spikein,
         total_reads_16s, total_reads_its, imte_reads, allo_reads, yarrowia_reads) %>%
  mutate(obs_allo = allo_reads / total_reads_16s,
         obs_imte = imte_reads / total_reads_16s,
         obs_its = yarrowia_reads / total_reads_its) %>%
  drop_na()
saveRDS(obs, "data/intermediate/obs.rds")
###2. expected proportion of spikein###
# 2.1. read qpcr data
qpcr <-
  readRDS("data/intermediate/qpcr.rds") %>%
  select(-sample) %>%
  group_by(marker, code) %>%
  dplyr::summarize_all(mean)

# 2.2. join qpcr and obs reads tables
# because it is needed to add to the qcpr the data from the spikein added
obs_vs_exp <-
  left_join(qpcr,
            obs,
            by = "code") %>%
  mutate(
    spikein = ifelse(marker == "16s", spikein * 5/7, # because the Kvf contained a mix of the its and 16s spikeins
                     ifelse(marker == "its", spikein * 2/7, NA)),
    kt_exp = kt_expected(So = s0, 
                         Kvf = ifelse(marker == "its", spikein * 10, # because I used a x10 lower dilution
                                      ifelse(marker == "16s", spikein, NA)),
                         Ko = k0, Kvo = kv0, Dk = dk, Ds = ds)) %>%
  select(marker, sample, code, s0, spikein, obs_allo:kt_exp) %>%
  pivot_longer(cols = c(obs_allo, obs_imte, obs_its),
               names_to = "spikein_id", values_to = "obs_prop") %>%
  filter(marker == "16s" & spikein_id == "obs_allo" |
           marker == "16s" & spikein_id == "obs_imte" |
           marker == "its" & spikein_id == "obs_its") %>%
  # in the spikein there are equal number of cells from imte and allo but their copy number varies
  # from 7 in allo and 3 in imte. So the kt_ext for each of the species varies according to their copies
  mutate(kt_exp = ifelse(spikein_id == "obs_allo", kt_exp * 7/10, # correct for 16s copies
                         ifelse(spikein_id == "obs_imte", kt_exp * 3/10, kt_exp)))# correct for 16s copies

saveRDS(obs_vs_exp, "data/intermediate/obs_vs_exp.rds")
