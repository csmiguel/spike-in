ps_filter_split_taxa_ps <- function(ps = NULL, asv2second_ps = NULL, nameSecond_ps = NULL) {
  #split ps according to vector of taxa names.
  # a second ps object is sent to the Global Environment named "nameSecond_ps" and cntaining ASVs "asv2second_ps".
  # it returns a phyloseq without "asv2second_ps" asvs.
  require(phyloseq)
  assign(x = nameSecond_ps,
         value = prune_taxa(asv2second_ps, ps),
         envir = .GlobalEnv)
  h <- prune_taxa(!taxa_names(ps) %in% asv2second_ps, ps)
  cat("\nTaxa",
      asv2second_ps, "were moved from the returned phyloseq",
      " to another phyloseq called",
      nameSecond_ps, ", sent to the Global Environment.\n")
  return(h)
}
