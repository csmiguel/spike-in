# compute alpha diversity on the matrix of relative abundaces
library(phyloseq)
library(tidyverse)

#read filtered phyloseq
ktu16s <-
  readRDS("data/intermediate/ps16s_abs.rds") %>%
  estimate_richness() %>%
  tibble::rownames_to_column("lib") %>%
  select(lib, Observed, Shannon)

ktuITS <-
  readRDS("data/intermediate/psITS_abs.rds") %>%
  estimate_richness() %>%
  tibble::rownames_to_column("lib") %>%
  select(lib, Observed, Shannon)

alphadiv <- list(ktu16s = ktu16s,
                 ktuITS = ktuITS) %>%
  do.call(what = "cbind") %>%
  select(-ktuITS.lib) %>%
  rename(lib = ktu16s.lib) %>%
  mutate_at(vars(contains("Shannon")), .funs = round, 2)

# save objects
saveRDS(alphadiv,
        "data/intermediate/alpha_div.rds")
