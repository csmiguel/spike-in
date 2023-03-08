library(phyloseq)
library(tidyverse)
library(patchwork)

# function for plotting
source("code/functions/stacked_div_barplots.r")

# read phyloseq 
ps_bact <-
  readRDS("data/intermediate/ps16s_abs.rds")
ps_fung <-
  readRDS("data/intermediate/psITS_abs.rds")

# plotting par
w = 9
h = 8
# plots
# rel abundances bact phylum
pbact_rel <-
  plot_rel_abundace(
  ps = ps_bact,
  ntop = 10,
  glomtax = "phylum",
  pal_col = "Paired",
  relative = F,
  bacteria = T
)
# ggsave("output/stacked_abs_16s_phylum.pdf",
#        pbact_rel,
#        width = w,
#        height = h)

# rel abundances bact phylum
pfung_rel <-
  plot_rel_abundace(
    ps = ps_fung,
    ntop = 10,
    glomtax = "class",
    pal_col = "Paired",
    relative = F,
    bacteria = F
  )
# ggsave("output/stacked_abs_its_class.pdf",
#        pfung_rel,
#        width = w,
#        height = h)

pall <-
  pbact_rel / pfung_rel +
  plot_annotation(tag_levels = 'A')
ggsave(filename = "output/stacked_alpha_div.pdf",
       pall,
       height = h,
       width = w)
