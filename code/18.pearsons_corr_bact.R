# pearson's correlation between spike-ins
library(tidyverse)
# read data
ds1 <-
  readRDS("data/intermediate/abs_abundances.rds") %>%
  drop_na
sink("output/pearsons_bact_spikein.txt")
cor.test(x = ds1$wild_16s_gsoilAllo, y = ds1$wild_16s_gsoilImt, method = "pearson")
sink()
