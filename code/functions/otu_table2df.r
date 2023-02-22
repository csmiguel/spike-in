otu_table2df <- function(ps = NULL) {
  require(phyloseq)
  # Extract abundance matrix from the phyloseq object
  otu1 <- as(otu_table(ps), "matrix")
  # transpose if necessary
  if(!taxa_are_rows(ps)){otu1 <- t(otu1)}
  # Coerce to data.frame
  as.data.frame(otu1)
}