library(tidyverse)
# Table 2
# sample, dnaextcode, raw16s, rawITS, filtered16s, filteredITS, wild16S, wildITS, yarrowia, allo, imte, absolutewild16scopies, absolutewildITS_rel_count, Richness (KTUs), Shannon H


# metadata
meta <- readRDS("data/intermediate/meta.rds")
# raw reads
raw_reads <- readRDS("data/intermediate/raw_reads.rds")
# summary reads
abs <-
   readRDS("data/intermediate/abs_abundances.rds") %>%
   rename(dna_ext = sample)
# [1] "sample"             "total_reads_16s"    "total_reads_its"    "allo_reads"         "imte_reads"        
# [6] "yarrowia_reads"     "wild_reads_16s"     "wild_reads_its"     "code"               "its"               
# [11] "ssu"                "conc_ng_ul"         "date"               "campaign"           "rep"               
# [16] "soil_mg"            "spikein"            "timepoint"          "allo_16s"           "imte_16s"          
# [21] "yarrowia_cel"       "wild_16s_gsoilAllo" "wild_16s_gsoilImt"  "wild_ITScell_gsoil" "wild_16s_gsoil_av" 
sample_codes <-
  readRDS("data/intermediate/samples.rds") %>%
  tibble::rownames_to_column("dna_ext") %>%
  select(code, sample)

# join
left_join(abs, sample_codes,
          by = "code") %>%
  left_join(raw_reads, by = c("dna_ext" = "lib")) %>%
  select(sample, dna_ext, # sample data
         raw_reads_16s, raw_reads_its, # raw reads
         total_reads_16s, total_reads_its, # filtered reads
         wild_reads_16s, imte_reads, allo_reads, # breakdown filtered bacteria
         wild_reads_its, yarrowia_reads, # breakdown filtered fungi
         
         )
absolutewild16scopies, absolutewildITS_rel_count, Richness (KTUs), Shannon H
