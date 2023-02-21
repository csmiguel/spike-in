# this function takes the transformed output tax table from idTaxa and removes instances of "unidentified"
clean_its_idTaxa <- function(idtaxa = NULL) {
  stopifnot("domain" == dimnames(idtaxa)[[2]][1])
  # add domain to tax table
  idtaxa[,"domain"] <- "Fungi"
  # remove assigned taxonomy from first instance of "unidentified" to lowest taxonomic rank
    apply(idtaxa, 1, function(asv) {
      first_instance_unidentified <- grep("unidenti", asv)[1]
      if(!is.na(first_instance_unidentified)) {
        asv[first_instance_unidentified:length(asv)] <- NA
      }
      asv
    }) %>% t()
    }
