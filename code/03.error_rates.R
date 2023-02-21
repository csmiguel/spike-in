# compute error rates with custom script for Novaseq
library(dplyr)
library(ShortRead)
library(dada2)

# smoothing function to estimate errors from NovaSeq data
source("code/functions/loesserrfun.r")

#load data
filt_path <- "data/intermediate/filtered"
n <- 10000 # number of reads to subsample to estimate the error rates.
sample_names_its <-
  dir(path = filt_path,
      pattern = "ITS.*R1") %>%
  gsub("ITS_", "", x = .) %>%
  gsub("_R.*$", "", x = .)
# sample names 16S
sample_names_16S <-
  dir(path = filt_path,
      pattern = "16S.*R1") %>%
  gsub("16S_", "", x = .) %>%
  gsub("_R.*$", "", x = .)
# list sample samples
sample_names <- list(ITS = sample_names_its,
                     "16S" = sample_names_16S)

#sampler for 16S and ITS. It subsamples fastq files
c("16S", "ITS") %>%
  lapply(function(locus) {
  seq_along(sample_names[[locus]]) %>%
    lapply(function(x) {
      samplei <- sample_names[[locus]][x] #sample name
      #path to files
      fqf <-
        dir(path = filt_path,
            pattern = paste0(locus, "_", samplei, "_R1"),
            full.names = TRUE)
      fqr <-
        dir(path = filt_path,
            pattern = paste0(locus, "_", samplei, "_R2"),
            full.names = TRUE)
      #forward
      #subsample reads
      f <- ShortRead::FastqSampler(fqf, n)
      subsamplef <-  ShortRead::yield(f)
      close(f)

      pathf <-
        file.path(filt_path,
                  paste0(locus, "_",
                         samplei,
                        "_F_filt_10000_fastq.gz"))
      ShortRead::writeFastq(subsamplef, pathf)
      #reverse
      r <- ShortRead::FastqSampler(fqr, n)
      subsampler <-  ShortRead::yield(r)
      close(r)

      pathr <-
        file.path(filt_path,
                  paste0(locus, "_",
                         samplei,
                         "_R_filt_10000_fastq.gz"))
      ShortRead::writeFastq(subsampler, pathr)
      })
    })

#lear error rates
# forward ITS
errF_its <-
  list.files(filt_path,
             pattern = "ITS.*F_filt_10000_fastq.gz",
             full.names = TRUE) %>%
  dada2::learnErrors(multithread = TRUE,
                     nbases = 1e8,
                     errorEstimationFunction = loessErrfun_mod4)
# reverse ITS
errR_its <-
  list.files(filt_path,
             pattern = "ITS.*R_filt_10000_fastq.gz",
             full.names = TRUE) %>%
  dada2::learnErrors(multithread = TRUE,
                     nbases = 1e8,
                     errorEstimationFunction = loessErrfun_mod4)
# forward 16S
errF_16s <-
  list.files(filt_path,
             pattern = "16S.*F_filt_10000_fastq.gz",
             full.names = TRUE) %>%
  dada2::learnErrors(multithread = TRUE,
                     nbases = 1e8,
                     errorEstimationFunction = loessErrfun_mod4)
# reverse 16S
errR_16s <-
  list.files(filt_path,
             pattern = "16S.*R_filt_10000_fastq.gz",
             full.names = TRUE) %>%
  dada2::learnErrors(multithread = TRUE,
                     nbases = 1e8,
                     errorEstimationFunction = loessErrfun_mod4)

#plot errors
plot_errF_16S <- dada2::plotErrors(errF_16s, nominalQ = TRUE)
plot_errR_16S <- dada2::plotErrors(errR_16s, nominalQ = TRUE)
plot_errF_its <- dada2::plotErrors(errF_its, nominalQ = TRUE)
plot_errR_its <- dada2::plotErrors(errR_its, nominalQ = TRUE)

ggplot2::ggsave("output/plot_errF_16S.pdf", plot_errF_16S)
ggplot2::ggsave("output/plot_errR_16S.pdf", plot_errR_16S)
ggplot2::ggsave("output/plot_errF_its.pdf", plot_errF_its)
ggplot2::ggsave("output/plot_errR_its.pdf", plot_errR_its)

#save errors
save(errF_16s, errR_16s,
     errF_its, errR_its,
     file = "data/intermediate/errors.Rdata")
