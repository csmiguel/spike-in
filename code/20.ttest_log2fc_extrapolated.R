#t-test to compare if there are differentes between the observed and predicted quantities
# of the different spikins
library(tidyverse)
# custom function foldchange
source("code/functions/log2foldchange.R")
# load data
qpcr_codes <-
  readRDS("data/intermediate/qpcr.rds")$code %>%
  unique()
# obs vs expected data
h <-
  readRDS("data/intermediate/obs_vs_exp.rds") %>%
  mutate(qpcr = # either from qpcr or extrapollated
           as.numeric(code %in% qpcr_codes) %>%
           plyr::mapvalues(
             from = c(0,1),
             to = c("extrapolated", "qpcr")
           ),
         log2fc = log2foldchange(obs_prop, kt_exp)) %>%
  select(marker, qpcr, spikein_id, log2fc)

sink("output/ttest_log2fc_pqcrVSextrapolated.txt")
# Allobacillus
cat("\nAllobacillus\n")
t.test(
  filter(h, spikein_id == "obs_allo" & qpcr == "qpcr")$log2fc,
  filter(h, spikein_id == "obs_allo" & qpcr == "extrapolated")$log2fc,
)
# Imtechella
cat("\nImtechella\n")
t.test(
  filter(h, spikein_id == "obs_imte" & qpcr == "qpcr")$log2fc,
  filter(h, spikein_id == "obs_imte" & qpcr == "extrapolated")$log2fc,
)
# Yarrowia
cat("\nYarrowia\n")
t.test(
  filter(h, spikein_id == "obs_its" & qpcr == "qpcr")$log2fc,
  filter(h, spikein_id == "obs_its" & qpcr == "extrapolated")$log2fc,
)
sink()

