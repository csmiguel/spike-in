#assign taxonomy using DECIPHER
# a more detailed description in http://benjjneb.github.io/dada2/tutorial.html
# DECIPHER::idTaxa uses a better performing
# classification algorithm https://doi.org/10.1186/s40168-018-0521-5.

assign_taxonomy_decipher <- function(seqtabdada2 = NULL,
                                     training_set = NULL,
                                     subset = NULL,
                                     nprocessors = 1) {
  #seqtabdada2, is a OTU like table produced by makeSequenceTable
  # (and possibly filtered, eg removeBimeraDenovo)
  #training_set, is an R object formatte for DECIPHER. Downloaded from:
  # "http://www2.decipher.codes/Classification/TrainingSets"
  #nprocessors, number of processors to use. In MAC the max is 1
  # see more info here: https://github.com/benjjneb/dada2/issues/1333
  # extract sequences from OTU object
  require(DECIPHER)
  require(dplyr)
  dna <-
    Biostrings::DNAStringSet(row.names(seqtabdada2))
  if(is.numeric(subset))
    dna <- dna[sample(seq_along(dna), subset)]
    # DNAStringSet from the ASVs
    ids <- DECIPHER::IdTaxa(dna,
                            training_set,
                            #for MACs, only 1 processor (see below)
                            processors = nprocessors,
                            strand = "top",
                            verbose = T)

  # ranks of interest
  ranks <- c("domain", "phylum", "class", "order",
             "family", "genus", "species")
  # Convert the output object of class "Taxa" to a matrix analogous to
  # the output from assignTaxonomy
  taxid <-
    ids %>%
    sapply(function(x) {
      m <- match(ranks, x$rank)
      taxa <- x$taxon[m]
      taxa[startsWith(taxa, "unclassified_")] <- NA
      taxa
    }) %>% t()
  colnames(taxid) <- ranks
  rownames(taxid) <- as.character(dna)
  taxid
  }
