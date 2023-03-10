library(tidyverse)

# table 1
t1 <-
  readRDS("data/intermediate/meta.rds") %>%
  drop_na %>%
  select(code, date, campaign, timepoint) %>%
  arrange(code, date) %>%
  distinct() %>%
  mutate(timepoint = stringr::str_replace(timepoint, "fin", "end"),
         sample = LETTERS[1:nrow(.)]) %>%
  select(sample, date, timepoint, campaign, code)
saveRDS(t1, "data/intermediate/t1.rds")

# t2
qpcr <-
  readRDS("data/intermediate/qpcr.rds") %>%
  select(code, marker, sample, sm, s0, kvf_g_soil) %>%
  as_tibble() %>%
  group_by(marker, sample) %>%
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

t2 <-
  left_join(qpcr$`16s`, qpcr$its,
            by = c("code", "sample", "sm"))
saveRDS(t2, "data/intermediate/t2.rds")
