ps_filter_contaminants <- function(ps = NULL, blanks = NULL, excl = NULL) {
  # ps is a phyloseq object
  # blank is a regex pattern matching blanks
  # excl, is a vector with samples to exclude from decontam calculations
  require(decontam)
  require(phyloseq)
  #assertions
  assertthat::assert_that(length(grep(blanks, sample_names(ps))) > 1)
  #blanks and internal controls
  blank_samples <- grep(blanks, sample_names(ps), value = T)
  blank_postions <- # logical vector
    sample_names(ps) %in% blank_samples %>%
    {.[rownames(otu_table(ps)) %in% excl] <- NA; .} # NA to internal controls
  #contanimanion analysis
  contam_df <- decontam::isContaminant(ps, neg = blank_postions)
  # contaminant ASVs
  contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])
  # non-contaminant ASVs
  non_contam_asvs <- row.names(contam_df[contam_df$contaminant == FALSE, ])
  h <-
    prune_taxa(non_contam_asvs, ps) %>%
    {prune_samples(
      !sample_names(.) %in% blank_samples, .)}
  if (identical(contam_asvs, character(0))) {
    cat("\nNo contaminant ASVs were identified.\n")
  } else if (length(contam_asvs) > 0) {
    cat("\n", contam_asvs,
        "was/were identified as contaminants and removed from the phyloseq object.\n")
  cat("\nSample/s",
      sample_names(h), "were kept.\nSample/s",
      blank_samples,
      "was/were removed.\n")
  }
  #report results from contaminant filtering
  no_reads <-
    mean(rowSums(otu_table(h)) /
           rowSums(otu_table(ps)[!sample_names(ps) %in% blank_samples]))
  no_taxa <- ntaxa(h) / ntaxa(ps)
  #report some results from the filtering
  cat("\nAfter removing contaminants and the blanks",
      blank_samples, ", a total of",
      no_reads, "of the reads and a",
      no_taxa, "of the ASVs were retained.\n")
  return(h)
}
