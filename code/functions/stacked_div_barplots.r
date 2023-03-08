plot_rel_abundace <- function(ps = NULL, ntop = NULL, glomtax = "phylum",
  pal_col = "Paired", relative = T, bacteria = T) {
  # ps, phyloseq object
  # ntop, number of taxa with highest abundances to plot
  # glomtax, rank to agglomerate by
  # pal_col, name of the palette color to use
  # relative, T/F: relative frequencies/absolute counts are plotted
  # bacteria, T/T: bacteria/fungi phyloseq as input
  # transform counts and agglomerate

  psx <-
    phyloseq::tax_glom(ps,
                       taxrank = glomtax,
                       NArm = FALSE)
  yylab <- "Bacteria 16S copies / g soil"
  if(!bacteria)
    yylab <- "Fungi units relative to spike-in / g soil"
  if(relative) {
    psx <-
      phyloseq::transform_sample_counts(psx, function(x) { # transform to proportions
        x / sum(x) })
    yylab <- "Relative abundance"
  }
  #
  
  if(phyloseq::taxa_are_rows(psx)) {
    h <- rowSums(otu_table(psx))
  } else if(!phyloseq::taxa_are_rows(psx)) {
    h <- colSums(otu_table(psx))
  }
  names_ntop <-
    h %>%
    sort(decreasing = T) %>%
    names() %>%
    .[1:ntop]
  #tax_table
  tax_table_rank <-
    tax_table(psx)@.Data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("asv")
  #vector with most abundant phyla
  names_legend <-
    c(tax_table_rank[[glomtax]][match(names_ntop, tax_table_rank$asv)],
      "other")
  # agglomerate low frequency taxa
  tidy_counts <- 
    phyloseq::psmelt(psx) %>%
    as_tibble() %>%
    select(Sample, OTU, Abundance, all_of(glomtax)) %>%
    mutate(labels_ntop = ifelse(OTU %in% names_ntop,
                                eval(parse(text = glomtax)),
                                "other")) %>%
    select(-all_of(glomtax)) %>%
    group_by(Sample, labels_ntop) %>%
    summarise(sum_tax = sum(Abundance)) %>%
    mutate(tax_col = factor(labels_ntop,
                            levels = names_legend)) %>%
    select(-labels_ntop)
  
  # relative abundance plots
  ggplot(tidy_counts) +
    geom_bar(aes(x = Sample, y = sum_tax, fill = tax_col),
             stat = "identity",
             position = "stack",
             color = "black",
             size = 0.1) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = -90, hjust = 0),
          strip.background = element_blank(),
          strip.text.x = element_text(size = 10, face = "italic")) +
    scale_fill_brewer(palette = pal_col, name = glomtax) +
    ylab(yylab)
  }
