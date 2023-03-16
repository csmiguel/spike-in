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
                 ktuITS = ktuITS)

# save objects
saveRDS(alphadiv,
        "data/intermediate/alpha_div.rds")
