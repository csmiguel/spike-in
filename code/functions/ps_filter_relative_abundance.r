ps_filter_relative_abundance <- function(ps = NULL,
                                      mean_prop = NULL, # recommended 1e-5
                                      n_samples = NULL) {
  # ps, phyloseq object.
  # mean_prop, average proportion of reads that a given ASV has across samples.
  # n_samples, number of samples which have to meet condition mean_prop.
  # if n_samples is NULL, ASVs are filtered taking into account only means.
  # if n_samples is numeric, ASVs are filtered with sum of mean_prop for at
  # least n_samples
  require(phyloseq)
  assertthat::assert_that(!phyloseq::taxa_are_rows(ps))
  #transform read count per sample to relative: sum == 1.
  ps_t <- transform_sample_counts(ps, function(x) x / sum(x))
  #filter ASV (get list of ASV to retain)
  if(is.null(n_samples)) {
    h <- filter_taxa(ps_t, function(x) mean(x) > mean_prop, TRUE)
  } else if (!is.null(n_samples)) {
    h <- filter_taxa(ps_t, function(x) sum(x > mean_prop) > n_samples, T)
        }
    #filter taxa
    hh <- prune_taxa(taxa_names(ps) %in% taxa_names(h), ps)
    #prop of reads kept
    no_reads <- mean(rowSums(otu_table(hh)) / rowSums(otu_table(ps)))
    no_taxa <- ntaxa(hh) / ntaxa(ps)    #report some results from the filtering
    cat("\nBy using", mean_prop, "as a threshold for mean proportion of reads",
        "across all samples, a proportion of",
        no_reads, "of the reads and a",
        no_taxa, "of the ASVs were retained.\n")
    return(hh)
    }
