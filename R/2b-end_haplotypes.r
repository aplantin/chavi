#' end_haplotypes
#'
#' Ends consensus risk haplotype for particular subjects (when fewer than
#' max.leave differ from majority)
#'
#' @keywords internal
#'
#' @param direction Which direction to move from site of interest. Options are "right" (genomic positions higher than site of interest) and "left" (genomic positions lower than site of interest).
#' @param consensus Data frame containing consensus haplotype segments.
#' @param riskhaps Vector of identifiers for risk haplotypes
#' @param haps Matrix of phased haplotypes
#' @param end.pos Current end position of consensus haplotype segment
#' @param ends Vector storing end.pos for all consensus haplotype segments
#' @param groups Internal variable tracking which subjects are in the current haplotype group
#' @param this.hapid Consensus haplotype segment currently under consideration
#' @return Consensus, groups, end.pos, ends
#' @export
#'
end_haplotypes <- function(direction, consensus, riskhaps, haps,
                           end.pos, ends, groups, this.hapid) {

  ## which are leaving?
  hap.tab <- table(as.numeric(haps[end.pos, groups[[this.hapid]]]))
  rare <- as.numeric(names(which.min(hap.tab)))
  leave.index <- which(haps[end.pos, groups[[this.hapid]]] == rare)
  leave.id <- colnames(haps[, groups[[this.hapid]]])[leave.index]

  which.ending <- consensus$id %in% leave.id & consensus$hapID == this.hapid
  if (direction == "left") {
    consensus[which.ending, ]$start <- haps$POS[end.pos]
  } else {
    consensus[which.ending, ]$end <- haps$POS[end.pos]
  }
  groups[[this.hapid]] <- groups[[this.hapid]][-leave.index]
  ends[[this.hapid]] <- end.pos

  return(list(consensus = consensus, groups = groups,
              end.pos = end.pos, ends = ends))
}

