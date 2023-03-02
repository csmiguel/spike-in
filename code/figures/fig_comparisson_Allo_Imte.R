# correlation between Allobacillus and Imtechella estimations
library(tidyverse)
library(ggpubr)

# read data
ds1 <-
  readRDS("data/intermediate/abs_abundances.rds") %>%
  drop_na
yrange <- range(ds1$wild_16s_gsoilAllo)
xrange <- range(ds1$wild_16s_gsoilImt)
all_lims <- range(c(xrange, yrange))

# correlation between 16s estimation using both bacteria
# add correlation estimates to plot
pcorr <-
  ggscatter(ds1, x = "wild_16s_gsoilImt", y = "wild_16s_gsoilAllo",
            add = "reg.line", conf.int = TRUE,
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "Inferred from Imtechella", ylab = "Inferred from Allobacillus") +
  stat_regline_equation(formula = y ~ 0 + x, label.y = 2*10^9) +
  geom_abline(slope = 1, intercept = 0, colour = "grey", linetype = "dashed") +
  xlim(all_lims* c(0.8, 1.1)) + ylim(all_lims* c(0.8, 1.1))
  
# save plots
ggsave(pcorr,
       filename = "output/corr_16s_spikes.pdf",
       width = 7, height = 6)
