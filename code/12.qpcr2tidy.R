# read qpcr data to tidy object
library(plyr)
library(tidyverse)
library(readxl)

source("code/functions/qpcr_template2tidy.R") # read sm template spikein excel
# sample codes from qpcr6
id_qpcr <- readRDS("data/intermediate/id_qpcr6.rds")

# read excel
filep <- "manuscript/sm_calcutaions_spikein.xlsx"
sheets <-
  readxl::excel_sheets(filep)
# read excel qpcr template into tidy format
qpcr <-
  lapply(sheets, function(x) {
    qpcr_template2tidy(filep, sheet = x) %>%
      select(-well) %>%
      mutate_at(vars(-"sample"), as.numeric) %>%
      mutate(sample = stringr::str_remove(sample, "_1/25")) %>%
      mutate(code = mapvalues( # assign code names to samples
        x = sample,
        from = id_qpcr$sample,
        to = id_qpcr$code),
        totals0 = s0 * ds / sm # total units of standard per gram of soil
      )
  }
  ) %>%
  setNames(sheets) %>%
  plyr::ldply(.id = "marker")

saveRDS(qpcr, "data/intermediate/qpcr.rds")
