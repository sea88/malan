#' Generate a function to simulate pedigree founder haplotype
#' 
#' @param ladder_min vector of minimum alleles; ladder_min[i] is the minimum allele at locus i
#' @param ladder_max vector of minimum alleles; ladder_max[i] is the maximum allele at locus i
#' @export
generate_get_founder_haplotype <- function(ladder_min, ladder_max) {
  stopifnot(is.vector(ladder_min))
  stopifnot(is.vector(ladder_max))
  stopifnot(length(ladder_min) == length(ladder_max))
  
  function() {
    unlist(lapply(seq_along(ladder_min), function(loc) {
      sample(seq(ladder_min[loc], ladder_max[loc]), 1)
    }))
  }
}

