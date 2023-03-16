library(ShortRead)
# correct from file name to sample name
correct_rownames <- function(x) {
  rn <- rownames(x) %>%
    gsub(pattern = "ITS_", replacement = "") %>%
    gsub(pattern = "16S_", replacement = "") %>%
    gsub(pattern = "_R.*$", replacement = "")
  rownames(x) <- rn
  return(x)
}
# number of raw reads
# path to filtered fastq
# 16S
raw16s <-
  list.files("data/raw",
             pattern = "16S.*R1.fastq.gz",
             full.names = TRUE) %>%
  ShortRead::countFastq() %>%
  dplyr::select(records) %>%
  dplyr::rename(raw_reads_16s = records) %>%
  correct_rownames() %>%
  tibble::rownames_to_column("lib")

rawITS <-
  list.files("data/raw",
             pattern = "ITS.*R1.fastq.gz",
             full.names = TRUE) %>%
  ShortRead::countFastq() %>%
  dplyr::select(records) %>%
  dplyr::rename(raw_reads_its = records) %>%
  correct_rownames() %>%
  tibble::rownames_to_column("lib")

raw <- full_join(raw16s, rawITS, by = "lib")

saveRDS(raw, "data/intermediate/raw_reads.rds")
