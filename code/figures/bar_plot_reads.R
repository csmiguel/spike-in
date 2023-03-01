# plot expected versus observed proportion of spikein
library(tidyverse)
obs <- readRDS("data/intermediate/abs_abundances.rds") %>%
  drop_na() %>%
  select(sample, code, wild_reads_16s, wild_reads_its, imte_reads, allo_reads, yarrowia_reads) %>%
  pivot_longer(cols = c(wild_reads_16s, wild_reads_its, imte_reads, allo_reads, yarrowia_reads),
              names_to = "nature", values_to = "reads") %>%
  mutate(nature = factor(nature,
                         levels = c("wild_reads_16s","wild_reads_its",
                                    "imte_reads", "allo_reads", "yarrowia_reads"))) %>%
  arrange(code, sample)
obs_bact <- obs %>%
  filter(nature %in% c("wild_reads_16s", "imte_reads", "allo_reads"))
obs_fung <- obs %>%
  filter(nature %in% c("wild_reads_its", "yarrowia_reads"))

ggplot() +
  geom_col(data = obs_bact,  aes(x = sample, y = -reads, fill = nature), position="stack") +
  geom_col(data = obs_fung,  aes(x = sample, y = reads, fill = nature), position="stack") +
  geom_vline(xintercept = 0) +
  scale_fill_manual(values = c("#2ca25f", "#b2e2e2", "#66c2a4" , "#fee8c8", "#fdbb84")) +
  coord_flip() +
  theme_classic() +
  theme(legend.title = element_blank())
