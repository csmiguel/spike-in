kt_expected <- function(So, Kvf, Ko, Kvo, Dk, Ds) {
  1- (So / (So + (Kvf * Ko * Dk / Kvo / Ds)))
}
