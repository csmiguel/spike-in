# kt expected verus observed
library(tidyverse)
# custom function foldchange
source("code/functions/log2foldchange.R")
myexp <- function(base, exp) {
  base^exp
}
# load data
obs_vs_exp <-
  readRDS("data/intermediate/obs_vs_exp.rds") %>%
  mutate(log2fc = log2foldchange(obs_prop, kt_exp)) %>%
  select(-c(ct:sm)) %>%
  group_by(spikein_id) %>%
  summarise(log2fc_av = mean(log2fc) %>% round(2),
            log2fc_sd = sd(log2fc)  %>% round(2)) %>%
  mutate(real = 2^log2fc_av  %>% round(2)) 
# write results
write.table(obs_vs_exp,
            "output/obs_vs_exted_log2foldchange.txt")
