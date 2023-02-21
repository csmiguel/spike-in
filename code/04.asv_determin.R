# asv determination
library(dplyr)
library(dada2)

# errors
load("data/intermediate/errors.Rdata")
source("code/parameters/filtering_ASVs.r")

# path to filtered fastq
filtpath <- "data/intermediate/filtered"
# 16S
filtF_16s <-
  sort(list.files(filtpath,
                  pattern = "16S.*R1.filt",
                  full.names = TRUE))
filtR_16s <-
  sort(list.files(filtpath,
                  pattern = "16S.*R2.filt",
                  full.names = TRUE))
# ITS
filtF_its <-
  sort(list.files(filtpath,
                  pattern = "ITS.*R1.filt",
                  full.names = TRUE))
filtR_its <-
  sort(list.files(filtpath,
                  pattern = "ITS.*R2.filt",
                  full.names = TRUE))

# ASV inference
dadaFs16S <-
  dada2::dada(filtF_16s,
             err = errF_16s,
             multithread = TRUE)
dadaRs16S <-
  dada2::dada(filtR_16s,
             err = errR_16s,
             multithread = TRUE)
dadaFsITS <-
  dada2::dada(filtF_its,
             err = errF_its,
             multithread = TRUE)
dadaRsITS <-
  dada2::dada(filtR_its,
             err = errR_its,
             multithread = TRUE)

# merge pairs
# 16S
merged_16S <-
  dada2::mergePairs(dadaFs16S, filtF_16s,
                   dadaRs16S, filtR_16s, verbose=TRUE)
names(merged_16S) <-
  sapply(strsplit(names(merged_16S), "_"), "[", 2)
# ITS
merged_ITS <-
  dada2::mergePairs(dadaFsITS, filtF_its,
                   dadaRsITS, filtR_its, verbose=TRUE)
names(merged_ITS) <-
  sapply(strsplit(names(merged_ITS), "_"), "[", 2)

#construct sequence table (analogous to an OTU table)
seqtab16S <- dada2::makeSequenceTable(merged_16S)
seqtabITS <- dada2::makeSequenceTable(merged_ITS)

#number of samples and variants
dim(seqtab16S)
dim(seqtabITS)
#distribution of sequence lengths
cat("Sequence length distribution of ASVs before removing chimeras")
cat("\nFor 16S:")
table(nchar(getSequences(seqtab16S)))
cat("\nWe retained reads in the range", range(filtering_seqtab_16S))
cat("\nFor ITS:")
cat("\nWe retained all reads since intronic region can be very variable in size")

#sequences that are much longer or shorter than expected may be the result of
# non-specific priming, and may be worth removing. I will only use this size filter for
# the 16S.
seqtab16S <-
  seqtab16S[, nchar(colnames(seqtab16S)) %in% filtering_seqtab_16S]

#remove quimeras
seqtabNoC_16S <- dada2::removeBimeraDenovo(seqtab16S)
seqtabNoC_ITS <- dada2::removeBimeraDenovo(seqtabITS)

#save objects
# dada
save(dadaFs16S, dadaRs16S,
     dadaFsITS, dadaRsITS,
     file = "data/intermediate/dada.Rdata")
#merged seqs
save(merged_16S, merged_ITS,
    file = "data/intermediate/mergers.Rdata")
# OTU like table from merged and 16S size-filtered sequences.
save(seqtab16S, seqtabITS,
    file = "data/intermediate/seqtab.Rdata")
#OTU tables no-chimeras
save(seqtabNoC_16S, seqtabNoC_ITS,
    file = "data/intermediate/seqtabNoC.Rdata")
