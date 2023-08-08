# filter ITS2 sequences meeting expected homology
library(tidyverse)
fp <- "data/intermediate/ITSx"
fp1 <- file.path(fp, "ITSx.positions.txt")

# output from ITSx summary
df1 <-
  read.delim(fp1, sep = "\t", skip = 0, header = F) %>%
  as_tibble() %>%
  setNames(c("Sequence ID", "Length of the sequence",
             "SSU range", "ITS1 range", "rrna5.8S range",
             "ITS2 range", "LSU range", "problems")) %>%
  janitor::clean_names() %>%
  mutate(problems = as.factor(problems)) %>%
  select(-ssu_range) %>%
  # filter sequences starting at 5.8S spanning over ITS2 and ending at LSU.
  filter(grepl("Not found", its1_range),# no seq on ITS1
         grepl("5.8S: No start", rrna5_8s_range), # starting at 5.8S
        !grepl(" 1-", its2_range), # not starting at ITS2
        !grepl("Not found", lsu_range)) # with LSU present

# read ITS2
fp3 <-
  file.path(fp, "ITSx.ITS2.fasta")
seqs_its2 <-
  Biostrings::readDNAStringSet(fp3)

# filter ITS2
seqs_its2_filt <-
  sapply(df1$sequence_id, function(x) grep(x, names(seqs_its2))) %>%
  {seqs_its2[.]}

# write ITS2 filtered
Biostrings::writeXStringSet(seqs_its2_filt, filepath = file.path(fp, "its2_itx.fasta"))

