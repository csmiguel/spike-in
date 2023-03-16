# create supplementary file with sequence information
library(tidyverse)
library(readxl)

# raw reads
raw_reads <-
  readRDS("data/intermediate/raw_reads.rds")
# alpha diversity
alpha_div <-
  readRDS("data/intermediate/alpha_div.rds")
# summary reads
abs <-
   readRDS("data/intermediate/abs_abundances.rds") %>%
   rename(dna_ext = sample)
# sample codes
sample_codes <-
  readRDS("data/intermediate/samples.rds") %>%
  select(code, sample)

# join data
all_data <-
  left_join(abs, sample_codes,
          by = "code") %>%
  left_join(raw_reads,
            by = c("dna_ext" = "lib")) %>%
  left_join(alpha_div,
            by = c("dna_ext" = "lib"))
# 16s data
bact_data <-
  all_data %>%
  select(sample, dna_ext, # sample data
       raw_reads_16s, # raw reads
       total_reads_16s, # filtered reads
       wild_reads_16s, imte_reads, allo_reads, # breakdown filtered bacteria
       wild_16s_gsoil_av, # absolute per gram of soil
       "ktu16s.Observed", "ktu16s.Shannon" # alpha div
       ) %>%
  drop_na %>%
  setNames(
    c("sample", "dna_ext", "raw_reads",
      "filtered_reads", "wild_reads",
      "Imtechella", "Allobacillus", "16s_gsample",
      "Richness", "Shannon")
  ) %>%
  as.data.frame()
# ITS data
its_data <-
  all_data %>%
  select(sample, dna_ext, # sample data
         raw_reads_its, # raw reads
         total_reads_its, # filtered reads
         wild_reads_its, yarrowia_reads, # breakdown filtered fungi
         wild_ITScell_gsoil, # absolute per gram of soil
         "ktuITS.Observed", "ktuITS.Shannon" # alpha div
  ) %>%
  drop_na %>%
  setNames(
    c("sample", "dna_ext", "raw_reads",
      "filtered_reads", "wild_reads","Yarrowia",
      "ITS_gsample",
      "Richness", "Shannon")
  ) %>%
  as.data.frame()
# variables explained
variab <-
  data.frame(
    variable = c("sample", "dna_ext", "raw_reads",
             "filtered_reads", "wild_reads",
             "Imtechella/Allobacillus/Yarrowia", "ITS/16s_gsample",
             "Richness", "Shannon"),
    description = c("sample name", "code for DNA extraction", "Raw sequencing reads",
             "Filtered reads", "Filtered reads corresponding to non-spike-in taxa",
             "Filtered reads assigned to any of the spike-in taxa",
             "Absolute 16S copies or ITS cell equivalents of wild taxa in 1 gram of the soil sample",
             "Alpha diversity: number of different KTUs", "Shannon H")
  )
# write excel
xlsx::write.xlsx(variab, row.names = F,
                 "output/sm3.xlsx", sheetName = "variables")
xlsx::write.xlsx(bact_data, row.names = T,
                 "output/sm3.xlsx", sheetName = "16s", append = T)
xlsx::write.xlsx(its_data, row.names = F,
                 "output/sm3.xlsx", sheetName = "ITS", append = T)
# save dfs
saveRDS(list(bact = bact_data,
             its = its_data),
        "data/intermediate/sm3.rds")
