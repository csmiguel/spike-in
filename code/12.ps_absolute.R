# transform counts for absolute counts in phyloseq objects
library(phyloseq)
library(tidyverse)

# load filtered phyloseq with no blanks and no internal controls
ps16s <-
  readRDS("data/intermediate/ps16s_ktu_filt_rel.rds")
psITS <-
  readRDS("data/intermediate/psITS_ktu_filt_rel.rds")

ab_ab <-
  readRDS("data/intermediate/abs_abundances.rds") %>%
  dplyr::filter(sample != "mcs168")

# transform sample counts according to absolute abundances
# NO CORRECTION for taxon-specific copy number is applied

# 16S
# create vector with size factors and order as sample names
# the result after correction will be the number of 16 copies in 1g of soil
size_factors16s_gsoil <- function() {
  m <- match(sample_names(ps16s), ab_ab$sample)
  # ratio copies in soil sample respect to reads
  r <- ab_ab$wild_16s_gsoil_av[m] / ab_ab$wild_reads_16s[m]
  return(r)
}

# Divide the counts in sample n by size_factors[n] and save in matrix `mat`
mat16s <-
  sweep(otu_table(ps16s), 1 + taxa_are_rows(ps16s), size_factors16s_gsoil(), FUN = "*")
# Note, the result is a plain matrix and not an `otu_table` object

# Make a new phyloseq object with the transformed counts
ps16s_abs <- ps16s
otu_table(ps16s_abs) <- otu_table(mat16s, taxa_are_rows = taxa_are_rows(ps16s))

# ITS
# create vector with size factors and order as sample names
# wild_copies / wild_reads: the result after correction will be the number of ITS-yarrowia equivalent copies in 1g of soil.
size_factorsITS_gsoil <- function() {
  m <- match(sample_names(psITS), ab_ab$sample)
  # ratio copies in soil sample respect to reads
  r <- ab_ab$wild_ITScell_gsoil[m] / ab_ab$wild_reads_its[m]
  return(r)
}
# Divide the counts in sample n by size_factors[n] and save in matrix `mat`
matits <-
  sweep(otu_table(psITS), 1 + taxa_are_rows(psITS), size_factorsITS_gsoil(), FUN = "*")
# Note, the result is a plain matrix and not an `otu_table` object

# Make a new phyloseq object with the transformed counts
psITS_abs <- psITS
otu_table(psITS_abs) <- otu_table(matits, taxa_are_rows = taxa_are_rows(psITS))

# save phyloseq object with corrected absolute counts
saveRDS(ps16s_abs, "data/intermediate/ps16s_abs.rds")
saveRDS(psITS_abs, "data/intermediate/psITS_abs.rds")
