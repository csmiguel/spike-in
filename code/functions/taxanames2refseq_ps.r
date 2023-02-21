taxanames2refseq_ps <- function(ps) {
  #ps is phyloseq object with an otu_table, sample_data, and tax_table
  #input ps has seqs as taxa_names
  #ps taxa_names are renamed as ASVxx instead
  #returns ps with refseq slot with dna seqs
  require(Biostrings); require(phyloseq)
  dna <- Biostrings::DNAStringSet(taxa_names(ps))
  names(dna) <- phyloseq::taxa_names(ps)
  ps <- phyloseq::merge_phyloseq(ps, dna)
  taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
  ps
  }
