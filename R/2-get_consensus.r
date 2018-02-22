#' get_consensus
#'
#' Finds the consensus risk haplotype shared between case subjects
#'
#' @param haps Matrix of phased haplotypes. Rows are markers, columns are copies of the chromosome (2 columns per subject). Must have a column called POS containing the genetic positions of the markers. Other column names should have the format SUBJ1.1 and SUBJ1.2 for the two haplotypes of SUBJ1.
#' @param riskhaps Vector of risk haplotype identifiers (one per case subject).
#' @param site.start Genomic position at start of site of interest.
#' @param site.end Genomic position at end of site of interest.
#' @param verbose Logical indicating whether progress updates should be printed. Default FALSE.
#' @return List containing data frame of consensus haplotype segments and vector of split positions.
#' @export
#'
get_consensus <- function(haps, riskhaps, site.start, site.end, verbose = F) {

  ## set up output file
  consensus <- data.frame(id = riskhaps)
  consensus$end = consensus$start = NA
  consensus$hapID <- 1

  ## set up other files
  groups <- list()
  groups[[1]] <- riskhaps

  ## this constant controls number of haplotypes that may leave
  ## without starting a new haplotype group
  MAX.LEAVE <- 1

  res1 <- consensus_inner(direction = "right", consensus, haps, riskhaps,
                          site.start, site.end, groups, MAX.LEAVE, verbose)
  res2 <- consensus_inner(direction = "left", consensus, haps, riskhaps,
                          site.start, site.end, groups, MAX.LEAVE, verbose)

  ## combine results
  out1 <- res1$consensus
  out1[1:length(riskhaps), ]$start <- res2$consensus[1:length(riskhaps), ]$start
  out2 <- res2$consensus[-c(1:length(riskhaps)), ]
  out2$hapID <- out2$hapID + max(res1$consensus$hapID) - 1

  ## return combined results
  out.consensus <- rbind(out1, out2)
  out.splits <- sort(c(res1$splits, res2$splits))

  return(list(consensus = out.consensus, splits = out.splits))
}






