# return feature table from KTU:klustering with clustered OTUs and with DNA sequences as rownames
klustering2otutable <- function(kobject = NULL) {
  otut <- kobject$KTU.table
  seqshash <- kobject$ReqSeq
  stopifnot(all(names(seqshash) == rownames(otut)))
  rownames(otut) <- seqshash
  return(otut)
}



