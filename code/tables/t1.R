library(tidyverse)
library(flextable)
library(officer)

# table 1
t1 <-
  readRDS("data/intermediate/samples.rds")
# format t1_2_word
ft_t1 <- flextable(t1)
# write ft to word file
flextable::save_as_docx(ft_t1,
                        path = "output/t1.docx")

