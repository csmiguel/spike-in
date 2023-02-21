read_length <- 250 # read length from sequencing run.
amplicon_ITS <- 350 # bp of amplicon in metabarcoding.
amplicon_16S <- 300 # bp of amplicon in metabarcoding.
# reads were 250nt, but after trimming primers they were around 230nt.
# expected targetted ITS amplicon is 350nt, although it can change considerable due to introns.
# expected targetted 16S amplicon is around 300nt.
# quality plots from binned NovaSeq qualities show a good overall quality at the 3' end.
# anyways, since we still have a lot of overlapping I will truncate reads to 220.
trunc_f <- 220 # truncation length for forward reads.
trunc_r <- 220 # truncation length for reverse reads.
expected_errors <- c(3, 3) #maximum expected_errors for forward and reverse reads
