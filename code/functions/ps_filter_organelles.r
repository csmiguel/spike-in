ps_filter_organelles <- function(ps = NULL) {
  cloro <- phyloseq::subset_taxa(ps, order %in% "Chloroplast")
  mito <- phyloseq::subset_taxa(ps, family %in% "Mitochondria")
  h <- phyloseq::subset_taxa(
    ps, !order %in% "Chloroplast" & !family %in% "Mitochondria")
  cat("\nA total of", ntaxa(cloro), "ASVs from Chloroplast and", ntaxa(mito),
    "ASVs from Mitochondria were removed from the phyloseq object\n")
  assertthat::assert_that(
  sum(
    grepl(
      x = tax_table(h),
      pattern = "[Cc]hloroplast|[Mm]itochondria|[Ee]ukary")
    ) == 0,
  msg = "Not all organelle DNA has been filtered.")
  cat("\nAny taxonomic rank matching of",
      "[Cc]hloroplast|[Mm]itochondria|[Ee]ukary has been removed\n")
  no_reads <- mean(rowSums(otu_table(h)) / rowSums(otu_table(ps)))
  no_taxa <- ntaxa(h) / ntaxa(ps)
  #report some results from the filtering
  cat("\nAfter removing ASVs from organelles, a",
  no_reads, "of the reads and a",
  no_taxa, "of the ASVs were retained.\n")
  return(h)
  }
