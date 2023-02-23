# eq 7bis

wild_copies_gsoil <- function(sample_wild_reads, spikein_reads, spikein_copies, soilg) {
  # Sr, sample_wild_reads, non-spike in taxa in reads from the phyloseq object
  # Kr, spikein_reads, spike in reads in the phyloseq object
  # Kc, spikein_copies, number of spikein copies added to the soil sample at DNA extraction
  # Sm, soilg, grams of soil in DNA extraction
  sample_wild_copies = sample_wild_reads / spikein_reads * spikein_copies
  # for 1 gram of soil
  per_g_soil <- sample_wild_copies / soilg
  return(per_g_soil)
}

##functions###
# function to compute c2: c1.v1 = c2.v2
c2 <- function(c1, v1, v2) {
  h <- c1*v1/v2
  return(h)
}
