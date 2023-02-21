ps_filter_prevalence <- function(ps = NULL,
                              mult_threshold = 0,
                              prevalence_threshold = 0,
                              abundance_threshold = 0) {
  #ps, is a phyloseq object
  #mult_threshold, is the threshold to remove low frequent variants
  #mult is vector originating from the product of reads supporting asv across
  # all samples. ej asv1: 23, 5, 50 , 0 = 23*5*50 = 5,750
  #prevalence_threshold, present in more than n samples
  #abundance_threshold, with at least n reads supporting the asv.
  require(phyloseq)
  require(dplyr)
  assertthat::assert_that(!phyloseq::taxa_are_rows(ps))
  h <-
    data.frame(prevalence = apply(otu_table(ps),
                                  ifelse(taxa_are_rows(ps), yes = 1, no = 2),
                                  function(x){sum(x > 0)}),
              abundance = taxa_sums(ps),
              #vector with product of all non-0 abundances
              mult = apply(otu_table(ps),
                           ifelse(taxa_are_rows(ps), yes = 1, no = 2),
                           function(x){prod(x[which(x > 0)])}),
              tax_table(ps)) %>%
    dplyr::filter(prevalence > prevalence_threshold,
                  abundance > abundance_threshold,
                  mult > mult_threshold) %>%
    rownames() %>%
      {prune_taxa(taxa = ., x = ps)}
  #prop of reads kept
  filt_asvs <- mean(rowSums(otu_table(h)) / rowSums(otu_table(ps)))
  #report some results from the filtering
  no_reads <- mean(rowSums(otu_table(h)) / rowSums(otu_table(ps)))
  no_taxa <- ntaxa(h) / ntaxa(ps)
  cat("\nBy using",
      mult_threshold, "as a threshold for product of abundances,",
      prevalence_threshold, "as a threshold for prevalence and",
      abundance_threshold, "as a threshold for abundace,",
      no_reads, "of the reads and a",
      no_taxa, "of the ASVs were retained.\n")
  return(h)
  }
