# plot expected versus observed proportion of spikein
library(tidyverse)
# custom function foldchange
source("code/functions/log2foldchange.R")
# load data
obs_vs_exp <- readRDS("data/intermediate/obs_vs_exp.rds")
qpcr_codes <-
  readRDS("data/intermediate/qpcr.rds")$code %>%
  unique()
# distinguish between samples that were in the qpcr and those for which I estimated kt_expectec
# based on average s0
obs_vs_exp$qpcr <- as.numeric(obs_vs_exp$code %in% qpcr_codes)
# plot
p1 <-
 obs_vs_exp %>%
  ggplot() +
  geom_jitter(aes(x = spikein_id, y = log2foldchange(obs_prop, kt_exp),
                  color = as.character(qpcr),
              #    shape = as.character(qpcr)
              ),
              alpha = 0.7,
              width = 0.2,
              size = 4) +
  geom_abline(intercept = 0, slope = 0, linetype = 2, linewidth = 0.5, color = "grey") +
  ylab("Log2 fold change of observed/expected reads") +
  xlab("")

p2 <-
  p1 +
  scale_x_discrete(breaks=c("obs_allo", "obs_imte", "obs_its"),
                   labels=c("Allobacillus\nhalotolerans (16S)",
                            "Imtechella\nhalotolerans (16S)",
                            "Yarrowia\nlipolytica (ITS2)")) +
  scale_color_manual(values = c("grey70", "grey20"),
                     labels = c("extrapolated estimations","real estimations")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = c(.2, .8),
        legend.background = element_rect(linewidth = 0.5,
                                         linetype="solid",
                                         colour = "grey30"))
ggsave("output/log2foldchange_obs_vs_expected.pdf",
       width = 6.5,
       height = 5)
