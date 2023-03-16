# coefficient of variation between triplicates for absolute counts
library(tidyverse)
df <-
  readRDS("data/intermediate/sm3.rds")
samples_code <-
  readRDS("data/intermediate/samples.rds") %>%
  select(sample, code)

df$bact %>%
      left_join(samples_code,
                by = "sample") %>%
      mutate(code = stringr::str_remove(code, ".$")) %>%
      group_by(sample) %>% 
      summarise(
        obsav = mean(`16s_gsample`),
        sd = sd(`16s_gsample`),
        n = n(),
        coef_var = sd / obsav
      ) %>%
  drop_na

df$its %>%
  left_join(samples_code,
            by = "sample") %>%
  mutate(code = stringr::str_remove(code, ".$")) %>%
  group_by(sample) %>% 
  summarise(
    obsav = mean(`ITS_gsample`),
    sd = sd(`ITS_gsample`),
    n = n(),
    coef_var = sd / obsav
  ) %>%
  drop_na
