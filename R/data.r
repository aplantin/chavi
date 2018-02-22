#' Simulated haplotype matrix in format for input into get.consensus().
#'
#' @format A matrix with 50 rows and 21 columns:
#' \describe{
#'   \item{POS}{Genomic position of marker}
#'   \item{SUBJID1.1}{Data for subject 1, haplotype 1.}
#'   ...
#' }
"haps"



#' Haplotype identifiers for risk haplotypes (i.e., haplotypes carrying
#' the consensus risk haplotype at the site of interest).
#'
#' @format A vector of length 10, with one haplotype identifier per subject.
"riskhaps"
