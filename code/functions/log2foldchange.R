log2foldchange <- function(num, den) {
  fc <- num/den
  l2fc <- log2(fc)
  return(l2fc)
}
