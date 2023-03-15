library(tidyverse)

# qpcr data
qpcr <-
  readRDS("data/intermediate/qpcr.rds") %>%
  select(code, marker, sample, sm, s0, kvf_g_soil) %>%
  as_tibble() %>%
  rename(dna_ext = sample) %>%
  group_by(marker, dna_ext) %>%
  mutate(
    s0av = mean(s0),
    s0sd = sd(s0),
    kvfav = mean(kvf_g_soil),
    kvfsd = sd(kvf_g_soil),
  ) %>%
  ungroup() %>%
  select(-s0, -kvf_g_soil) %>%
  distinct() %>%
  plyr::dlply(~marker) %>%
  lapply(function(x) {
    names(x) <-
      paste0(
        names(x), 
        c(rep("",4), 
          rep(as.character(x$marker[1]), 4)
        )
      )
    select(x, -marker)
  })
# transform t1
samples_temp <-
  readRDS("data/intermediate/samples.rds") %>%
  tibble::rownames_to_column("dna_ext") %>%
  select(code, sample)

# format t2
t2 <-
  left_join(qpcr$`16s`, qpcr$its,
            by = c("code", "dna_ext", "sm")) %>%
  left_join(samples_temp,
            by = "code") %>%
  mutate(coef_var_16s = kvfsd16s / kvfav16s,
         coef_var_its = kvfsdits / kvfavits) %>%
  arrange(sample)
# coef of variations for technical replicas
range(t2$coef_var_16s)
range(t2$coef_var_its)

# summarize by template
t2_by_sample <-
  t2 %>%
  group_by(sample) %>%
  summarise(sample,
            kvfav_16s = mean(kvfav16s) %>% round(1),
            kvfsd_16s = sd(kvfav16s) %>% round(1),
            kvfav_its = mean(kvfavits / 10) %>% round(1),
            kvfsd_its = sd(kvfavits / 10) %>% round(1)
  ) %>%
  distinct()

coeff_trip16s <-
  t2_by_sample %>%
  mutate(coef_trip16s = kvfsd_16s / kvfav_16s,
         coef_tripits = kvfsd_its / kvfav_its)
# save
saveRDS(coeff_trip16s, "data/intermediate/coeff_qpcr_triplicates.rds")
