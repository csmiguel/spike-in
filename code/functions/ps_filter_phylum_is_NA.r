ps_filter_phylum_is_NA <- function(ps = NULL) {
  require(phyloseq)
 #filter out ASVs with no Phylum assigned
  h <- subset_taxa(ps, !is.na(phylum))
  no_reads <- mean(rowSums(otu_table(h)) / rowSums(otu_table(ps)))
  no_taxa <- ntaxa(h) / ntaxa(ps)
  #report some results from the filtering
  cat("\nAfter removing ASVs with no phylum assigned, a",
  no_reads, "of the reads and a",
  no_taxa, "of the ASVs were retained.\n")
  return(h)
  }
