#' pairwise_summary
#'
#' Summarizes overall start and end position of portion of
#' consensus risk haplotype shared between two particular
#' case individuals.
#'
#' @param consensus Data frame with all consensus haplotype segments for each case subject (output from get.consensus)
#' @param id1 First risk haplotype identifier (should have form like SUBJ1.2 and be an element of the riskhaps vector)
#' @param id2 Second risk haplotype identifier (should have form like SUBJ3.1 and be an element of the riskhaps vector)
#' @return Start, end, and length of haplotype shared between these subjects
#' @export
### Longest segment shared with a *single* other risk subject
pairwise_summary <- function(consensus, id1, id2) {
  this.i <- consensus[consensus$id == id1, ]
  this.j <- consensus[consensus$id == id2, ]
  shared.haps <- this.i$hapID[which(this.i$hapID %in% this.j$hapID)]

  ## on shared haps, find max end and min start
  this.i <- this.i[this.i$hapID %in% shared.haps, ]
  this.j <- this.j[this.j$hapID %in% shared.haps, ]
  pair.end <- min(max(this.i$end), max(this.j$end))
  pair.start <- max(min(this.i$start), min(this.j$start))

  ## start, end, length
  return(list(start = pair.start, end = pair.end,
              length = (pair.end - pair.start)))

}
