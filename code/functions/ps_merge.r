# set of functions to merge phyloseq objects. They are specific to this dataset
# merge otu tables with shared samples/asvs from 2 phyloseq objects
merge_otutables <- function(ps1 = NULL, ps2 = NULL) {
  # df with otu table from ps1. rows are samples
  ps11 <-
    otu_table2df(ps1) %>%
    t %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample")
  # df with otu table from ps2. rows are samples
  ps22 <-
    otu_table2df(ps2) %>%
    t %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample")
  # join otu tables
  newotu <- dplyr::left_join(ps11, ps22, by = "sample")
  # introduce 0s where NAs.
  newotu[is.na(newotu)] <- 0
  tibble::column_to_rownames(newotu, "sample")
}
merge_taxtable <- function(ps1 = NULL, ps2 = NULL, ps3 = NULL) {
  tax <- c(taxa_names(ps1), taxa_names(ps2))
  h <- prune_taxa(tax, ps3)
  tax_table(h)
}
