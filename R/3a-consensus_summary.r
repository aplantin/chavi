#' consensus_summary
#'
#' Summarizes overall start and end position of consensus risk haplotype
#' for each case individual (on their risk haplotype).
#'
#' @param consensus Data frame with all consensus haplotype segments for each case subject (output from get.consensus)
#' @return Data frame containing risk haplotype identifier and overall start/end of consensus haplotype
#' @export
consensus_summary <- function(consensus) {
  out.min <- aggregate(x = consensus$start, by = list(consensus$id), FUN = min)
  out.max <- aggregate(x = consensus$end, by = list(consensus$id), FUN = max)
  out.tab <- merge(out.min, out.max, by = "Group.1")
  colnames(out.tab) <- c("id", "start", "end")
  return(out.tab)
}

