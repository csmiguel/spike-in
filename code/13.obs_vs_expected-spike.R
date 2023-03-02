# plot expected versus observed proportion of spikein
library(tidyverse)
source("code/functions/kt_expected.R")

###1. observed proportion of spikein###
# proportion of observed spikein reads respect to total
obs <-
  readRDS("data/intermediate/abs_abundances.rds") %>%
  select(sample, code, spikein, soil_mg,
         total_reads_16s, total_reads_its, imte_reads,
         allo_reads, yarrowia_reads) %>%
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
  left_join(obs,qpcr,
            by = "code") %>%
  #select cases which were not in qpcr
  {
    df1 <- .[!(.$code %in% qpcr$code), ];
    # duplicate cases which are not in qpcr
    rbind(df1, .)
  } %>% arrange(marker, sample) %>%
  {  # replace NAs in other variables with values
    lnas <-sum(is.na(.$marker))
    .$marker[is.na(.$marker)] <- rep(c("16s", "its"), lnas/2)
    .$kt[is.na(.$kt)] <- rep(c(0.02, 0.01), lnas/2)
    .$kv0[is.na(.$kv0)] <- rep(c(10, 20), lnas/2)
    .$ds[is.na(.$ds)] <- .$ds[1]
    .$dk[is.na(.$dk)] <- .$dk[1]
    .$k0[is.na(.$k0)] <- .$k0[1]
    .
    
  } %>%
  mutate( # impute s0 data for NAs
    totals0 = ifelse(
      is.na(s0) & marker == "16s",
      mean(
        pull(filter(., marker == "16s" & !is.na(s0)), totals0)
      ),
      ifelse(
        is.na(s0) & marker == "its",
        mean(
          pull(filter(., marker == "its" & !is.na(s0)), totals0)
        )
        ,totals0)),
    s0 = ifelse(is.na(s0), totals0 * soil_mg /ds / 1000, s0),
    # because the Kvf contained a mix of the its and 16s spikeins
    spikein = ifelse(marker == "16s", spikein * 5 / 7,
                     ifelse(marker == "its", spikein * 2 / 7, NA)),
    kt_exp = kt_expected(So = s0,
                         # because I used a x10 lower dilution
                         Kvf = ifelse(marker == "its", spikein * 10,
                                      ifelse(marker == "16s", spikein, NA)),
                         Ko = k0, Kvo = kv0, Dk = dk, Ds = ds)) %>%
  select(marker, sample, code, s0, spikein, obs_allo:kt_exp) %>%
  pivot_longer(cols = c(obs_allo, obs_imte, obs_its),
               names_to = "spikein_id", values_to = "obs_prop") %>%
  filter(marker == "16s" & spikein_id == "obs_allo" |
           marker == "16s" & spikein_id == "obs_imte" |
           marker == "its" & spikein_id == "obs_its") %>%
  # in the spikein there are equal number of cells from imte and allo but
  # their copy number varies from 7 in allo and 3 in imte. So the kt_ext for
  # each of the species varies according to their copies
  # correct for 16s copies
  mutate(kt_exp = ifelse(spikein_id == "obs_allo", kt_exp * 7 / 10,
                         # correct for 16s copies
                         ifelse(spikein_id == "obs_imte", kt_exp  * 3 / 10, kt_exp)))


saveRDS(obs_vs_exp, "data/intermediate/obs_vs_exp.rds")
