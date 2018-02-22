#' split_haplotypes
#'
#' Splits the current consensus risk haplotype segment into two new groups
#'
#' @keywords internal
#'
#' @param direction Which direction to move from site of interest. Options are "right" (genomic positions higher than site of interest) and "left" (genomic positions lower than site of interest).
#' @param consensus Skeleton data frame to which consensus haplotype segments will be added
#' @param haps Matrix of phased haplotypes
#' @param groups Internal variable tracking which subjects are in the current haplotype group
#' @param this.hapid Consensus haplotype segment currently under consideration
#' @param end.pos Current end position of consensus haplotype segment
#' @param splits Vector of positions at which consensus haplotype splits into two new segments
#' @param ends Vector storing end.pos for all consensus haplotype segments
#' @return Consensus, groups, splits, ends
#' @export
#'
#'
split_haplotypes <- function(direction = c("right", "left"), consensus, haps,
                             groups, this.hapid, end.pos, splits, ends) {

  splits <- switch(direction,
                   right = c(splits, haps$POS[end.pos]),
                   left = c(splits, haps$POS[end.pos]))
  which.ending <- which(consensus$id %in% groups[[this.hapid]] &
                          consensus$hapID == this.hapid)
  if (direction == "right") {
    consensus[which.ending, ]$end <- haps$POS[end.pos]
  } else {
    consensus[which.ending, ]$start <- haps$POS[end.pos]
  }

  g0 <- which(haps[end.pos, groups[[this.hapid]]] == 0)
  g0.hid <- length(groups) + 1
  g1 <- which(haps[end.pos, groups[[this.hapid]]] == 1)
  g1.hid <- length(groups) + 2
  ends[[g0.hid]] = ends[[g1.hid]] = end.pos

  groups[[g0.hid]] <- colnames(haps[, groups[[this.hapid]]])[g0]
  groups[[g1.hid]] <- colnames(haps[, groups[[this.hapid]]])[g1]

  new.rows <- switch(direction,
                     right = data.frame(id = groups[[this.hapid]],
                                        start = haps$POS[end.pos],
                                        end = NA, hapID = NA),
                     left = data.frame(id = groups[[this.hapid]],
                                       start = NA, end = haps$POS[end.pos],
                                       hapID = NA))
  new.rows$hapID[new.rows$id %in% groups[[g0.hid]]] <- g0.hid
  new.rows$hapID[new.rows$id %in% groups[[g1.hid]]] <- g1.hid

  consensus <- rbind(consensus, new.rows)
  return(list(consensus = consensus, groups = groups,
              splits = splits, ends = ends))
}

