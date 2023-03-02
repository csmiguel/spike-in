# plot expected versus observed proportion of spikein
library(tidyverse)
source("code/functions/log2foldchange.R")

obs_vs_exp <- readRDS("data/intermediate/obs_vs_exp.rds")

p1 <-
 obs_vs_exp %>%
  #filter(marker == "16s") %>%
  ggplot() +
  geom_jitter(aes(x = spikein_id, y = log2foldchange(obs_prop, kt_exp), color = code), width = 0.2) +
  geom_abline(intercept = 0, slope = 0, linetype = 2, linewidth = 0.5, color = "grey") +
  ylab("Log2 fold change of observed/expected reads") +
  xlab("")

p2 <-
  p1 +
  scale_x_discrete(breaks=c("obs_allo", "obs_imte", "obs_its"),
                   labels=c("Allobacillus\nhalotolerans (16S)",
                            "Imtechella\nhalotolerans (16S)",
                            "Yarrowia\nlipolytica (ITS2)")) +
  theme_classic()

ggsave("output/log2foldchange_obs_vs_expected.pdf",
       width = 6.5,
       height = 5)
