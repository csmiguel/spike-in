# construct phyloseq post KTU computation
library(phyloseq)
library(tidyverse)

# copy ASV DNA seq to refseq slot and simplify otu_table names
source("code/functions/taxanames2refseq_ps.r")
# ktu object (list) to input for idTaxa (matrix)
source("code/functions/klustering2otutable.r")

# load seqs with taxonomy
taxid_its <- readRDS("data/intermediate/ITSx/taxid_its2.rds")

# load OTU tables
ktuits <-
  readRDS("data/intermediate/ITSx/ktu_its2.rds") %>%
  klustering2otutable()

# meta
meta <- readRDS("data/intermediate/meta.rds")

# assert all samples in the OTU table are present in metadata
assertthat::assert_that(
  all(colnames(ktuits) %in% rownames(meta)))


ps_itsktu <-
  phyloseq::phyloseq(
    t(otu_table(ktuits, taxa_are_rows = TRUE)),
    sample_data(meta),
    tax_table(taxid_its)) %>%
  taxanames2refseq_ps()

#save phyloseq objects
saveRDS(ps_itsktu, file = "data/intermediate/ITSx/psKTU.rds")
