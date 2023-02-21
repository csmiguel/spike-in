# filter and truncate reads
library(dplyr)
library(dada2)

# path to trimmed sequences from cutadapt
path <- "data/intermediate"
#load filtering parameters
source("code/parameters/trim.r")
#create path for sequences that will be filtered downstream
filt_path <- file.path(path, "filtered")
if (!file_test("-d", filt_path)) dir.create(filt_path)
# name loci
loci <- c("16S", "ITS")

# truncated reads in and out. For ITS I will not truncate reads as I do not
# know the size of the target DNA.
truncated_in_out <-
  loci %>%
  lapply(function(x) {
  # path to fastq files
  fnFs <-
    sort(list.files(path,
                    pattern = paste0("*", x, ".*R1.cutadapt.fastq.gz"),
                    full.names = TRUE))
  fnRs <-
    sort(list.files(path,
                    pattern = paste0("*", x, ".*R2.cutadapt.fastq.gz"),
                    full.names = TRUE))
  #path to filtered files
  filtFs <-
    gsub("intermediate/",
       "intermediate/filtered/",
       fnFs) %>%
    gsub("cutadapt", "filt", x = .)
  filtRs <-
    gsub("intermediate/",
         "intermediate/filtered/",
         fnRs) %>%
    gsub("cutadapt", "filt", x = .)

  #truncate reads
  if(x == "16S") {
      dada2::filterAndTrim(
        fnFs, filtFs,
        fnRs, filtRs,
        maxN = 0,
        maxEE = expected_errors,
        truncQ = 2,
        minLen = 50,
        compress = TRUE,
        multithread = TRUE,
        truncLen = c(trunc_f, trunc_r))
  } else if(x == "ITS") {
      dada2::filterAndTrim(
        fnFs, filtFs,
        fnRs, filtRs,
        maxN = 0,
        maxEE = expected_errors,
        truncQ = 2,
        minLen = 50,
        compress = TRUE,
        multithread = TRUE)
      }
  }) %>%
  setNames(loci)

#save objects
saveRDS(truncated_in_out, "data/intermediate/truncated_in_out.rds")
