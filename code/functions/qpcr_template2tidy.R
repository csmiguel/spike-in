# read data from excel template in SMx as tidy table
qpcr_template2tidy <- function(path2excel, ...) {
  require(janitor)
  require(tidyverse)
  require(readxl)
  options(warn = -1)
  readxl::read_xlsx(filep, ...) %>%
    {
      .[grep("well", pull(., 1)):nrow(.), ]
    } %>%
    janitor::row_to_names(1) %>%
    janitor::clean_names(.) %>%
    {
      .[, 1:grep("kvf_g_soil", names(.))]
    } %>%
    drop_na()
  }

