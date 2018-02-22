#' end_pair
#'
#' Ends consensus risk haplotype for a pair of subjects when they are the
#' last remaining pair in a consensus haplotype segment
#'
#' @keywords internal
#'
#' @param direction Which direction to move from site of interest. Options are "right" (genomic positions higher than site of interest) and "left" (genomic positions lower than site of interest).
#' @param consensus Skeleton data frame to which consensus haplotype segments will be added
#' @param haps Matrix of phased haplotypes
#' @param this.hapid Consensus haplotype segment currently under consideration
#' @param groups Internal variable tracking which subjects are in the current haplotype group
#' @param ends Vector storing end.pos for all consensus haplotype segments
#' @return Consensus, groups, ends
#' @export
#'
end_pair <- function(direction = c("right", "left"), consensus, haps,
                     this.hapid, groups, ends) {

  ## end at next difference because fewer than max.leave left
  ## (should only apply if < max.leave haps ended and there are only 2 left)

  end.pos <- ends[[this.hapid]]
  which.in <- which(colnames(haps) %in% groups[[this.hapid]])

  ## allele differences within group, current position to end of chromosome
  pos.remaining <- switch(direction,
                          right = haps[-c(1:end.pos), which.in],
                          left = haps[c(1:end.pos), which.in])
  this.same <- apply(pos.remaining, 1, FUN = function(x) length(unique(as.numeric(x))) == 1)

  if (all(this.same)) {
    end.pos <- switch(direction,
                      right = length(haps$POS),
                      left = 1)
  } else {
    end.pos <- switch(direction,
                      right = end.pos + min(which(this.same == FALSE)),
                      left = max(which(this.same == FALSE)))
  }
  leave.id <- colnames(haps[, groups[[this.hapid]]])
  which.ending <- consensus$id %in% leave.id & consensus$hapID == this.hapid
  if (direction == "right") {
    consensus[which.ending, ]$end <- haps$POS[end.pos]
  } else {
    consensus[which.ending, ]$start <- haps$POS[end.pos]
  }
  groups[[this.hapid]] <- groups[[this.hapid]][-c(1:2)]
  ends[[this.hapid]] <- end.pos

  return(list(consensus = consensus, groups = groups, ends = ends))
}
