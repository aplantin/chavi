#' consensus_inner
#'
#' Finds the consensus risk haplotype shared between case subjects
#' (main body of algorithm for get.consensus function).
#'
#' @keywords internal
#'
#' @param direction Which direction to move from site of interest. Options are "right" (genomic positions higher than site of interest) and "left" (genomic positions lower than site of interest).
#' @param consensus Skeleton data frame to which consensus haplotype segments will be added
#' @param haps Matrix of phased haplotypes
#' @param riskhaps Vector of identifiers for risk haplotypes
#' @param site.start Genomic position at start of site of interest
#' @param site.end Genomic position at end of site of interest
#' @param groups Internal variable tracking which subjects are in the current haplotype group
#' @param max.leave Maximum number of haplotypes that may leave a haplotype group without starting their own new consensus haplotype segment
#' @param verbose Logical indicating whether progress updates should be printed. Default FALSE.
#' @return List containing data frame of consensus haplotype segments and vector of split positions.
#' @export
#'
consensus_inner <- function(direction = c("right", "left"), consensus, haps,
                            riskhaps, site.start, site.end, groups,
                            max.leave, verbose) {

  # Constants
  npos = length(haps$POS)

  # Output set-up
  splits <- c()

  # Internal variables
  ends <- list()
  ends[[1]] <- switch(direction,
                      right = max(which(haps$POS <= site.end)),
                      left = min(which(haps$POS >= site.start)))
  this.hapid <- 1
  if (verbose) {
    print(paste("Hap", this.hapid, "to the", direction, "out of",
                length(groups), "total."))
  }

  # Find consensus risk haplotype
  while (this.hapid <= length(groups)) {
    end.pos <- ends[[this.hapid]]
    end.this.hap <- FALSE

    while (length(groups[[this.hapid]]) > 2 & !end.this.hap) {
      # indices of haplotypes in this group
      which.in <- which(colnames(haps) %in% groups[[this.hapid]])

      # this is TRUE if all SNPs at the site are same
      pos.to.check <- switch(direction,
                             right = haps[-c(1:end.pos), which.in],
                             left = haps[c(1:end.pos), which.in])
      this.same <- apply(pos.to.check, 1,
                         FUN = function(x) {
                           length(unique(as.numeric(x))) == 1
                         })

      if (all(this.same)) {
        ## No differences until end of chromosome
        leave.id <- colnames(haps[, groups[[this.hapid]]])
        which.leaving <- which(consensus$id %in% leave.id &
                                 consensus$hapID == this.hapid)
        if (direction == "right") {
          consensus[which.leaving, ]$end <- haps$POS[length(haps$POS)]
        } else {
          consensus[which.leaving, ]$start <- haps$POS[1]
        }
        this.hapid <- this.hapid + 1
        if (verbose) {
          print(paste("Hap", this.hapid, "to the", direction,
                      "out of", length(groups), "total."))}
      } else {
        ## At least one SNP differs between haplotypes before end of chromosome
        # find first site where difference exists
        end.pos <- switch(direction,
                          right = end.pos + min(which(!this.same)),
                          left = max(which(!this.same)))

        num.differ <- min(table(as.numeric(haps[end.pos, groups[[this.hapid]]])))
        if (num.differ <= max.leave) {
          # some haplotypes END because difference is in fewer than MAX.LEAVE
          temp1 <- end_haplotypes(direction, consensus, riskhaps, haps,
                                  end.pos, ends, groups, this.hapid)
          consensus <- temp1$consensus
          groups <- temp1$groups
          ends <- temp1$ends
          end.pos <- temp1$end.pos
        } else {
          temp2 <- split_haplotypes(direction, consensus, haps, groups,
                                    this.hapid, end.pos, splits, ends)
          consensus <- temp2$consensus
          groups <- temp2$groups
          splits <- temp2$splits
          ends <- temp2$ends
          this.hapid <- this.hapid + 1
          end.this.hap <- TRUE
          if (verbose == TRUE & this.hapid <= length(groups)) {
            print(paste("Hap", this.hapid, "to the", direction,
                        "out of", length(groups), "total."))
          }
        }

      } # ends else mismatch before end of chromosome
    } # ends while(length(groups[[this.hapid]]) > 2)

    if (length(groups[[this.hapid]]) == 2 & !end.this.hap) {
      # end haplotype where these two haplotypes differ
      temp3 <- end_pair(direction, consensus, haps, this.hapid, groups, ends)
      consensus <- temp3$consensus
      groups <- temp3$groups
      ends <- temp3$ends
      this.hapid <- this.hapid + 1
      if (verbose == TRUE & this.hapid <= length(groups)) {
        print(paste("Hap", this.hapid, "to the", direction,
                    "out of", length(groups), "total."))
      }
    }
  } # ends while(this.hapid <= length(groups))

  return(list(consensus = consensus, splits = splits))
}


