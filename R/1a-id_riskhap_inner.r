#' id_riskhap_inner
#'
#' Identify which haplotype is the risk haplotype (at the site of interest)
#' for each case individual. Called by id.riskhap().
#'
#' @keywords internal
#'
#' @param ibs IBS data.frame
#' @param case.ids Vector of case identifiers
#' @param site.start Genomic position at start of site of interest
#' @param site.end Genomic position at end of site of interest
#' @return A vector of risk haplotype IDs
#' @export
#'
id_riskhap_inner <- function(ibs, case.ids, site.start, site.end) {

  ## find IBS segments that are shared between cases and include site of interest
  ibs <- ibs[ which(ibs$id1 %in% case.ids & ibs$id2 %in% case.ids), ]
  site.ibs <- ibs[ which(ibs$start < site.start & ibs$end > site.end), ]

  ## list containing, for each case haplotype, a data frame with length of each
  ## IBS segment this hap shares with each other case hap at site of interest
  segs.hap1 = segs.hap2 = list()
  for (i in 1:length(case.ids)) {
    segs.hap1[[i]] = segs.hap2[[i]] = data.frame(hap1 = rep(0, length(case.ids)),
                                                 hap2 = rep(0, length(case.ids)))
  }
  sorted.segs = list(segs.hap1, segs.hap2)

  ## fill it in from site.ibs (each IBS segment covering mutation is categorized)
  for (i in 1:nrow(site.ibs)) {
    this.seg <- site.ibs[i, ]
    case1.index <- which(case.ids == this.seg$id1)
    case2.index <- which(case.ids == this.seg$id2)
    h1 <- this.seg$hap1
    h2 <- this.seg$hap2

    ## subject 1 (ID and haplotype) defines which data frame we're in
    ## subject 2 defines row and column of that data frame
    sorted.segs[[h1]][[case1.index]][case2.index, h2] <- this.seg$length

    ## now the converse
    sorted.segs[[h2]][[case2.index]][case1.index, h1] <- this.seg$length
  }

  ## For each case haplotype, find longer segment shared with each partner
  ##     (of the two partner haplotypes)
  ## Results in a nested list (level 1 = case haplotype, level 2 = case ID,
  ##     element = vector of max IBS with each other case subject at site)
  for (i in 1:length(case.ids)) {
    sorted.segs[[1]][[i]] <- apply(sorted.segs[[1]][[i]][-i, ], 1,
                                   FUN = function(partner) pmax(partner[1], partner[2]))
    sorted.segs[[2]][[i]] <- apply(sorted.segs[[2]][[i]][-i, ], 1,
                                   FUN = function(partner) pmax(partner[1], partner[2]))
  }

  ## vector of median length for each case haplotype
  med.h1 <- unlist(lapply(sorted.segs[[1]], FUN = function(x) median(x)))
  med.h2 <- unlist(lapply(sorted.segs[[2]], FUN = function(x) median(x)))
  rh <- as.numeric(med.h1 < med.h2) + 1

  ## use mean as tiebreaker
  mean.h1 <- unlist(lapply(sorted.segs[[1]], mean))
  mean.h2 <- unlist(lapply(sorted.segs[[2]], mean))
  if(any(med.h1 == med.h2)) {
    eqs <- which(med.h1 == med.h2)
    warning(paste("Medians are equal for ",
                  paste(case.ids[eqs], collapse = " and "),
                  "; using mean instead", sep = ""))
    uneq <- c()

    for (i in 1:length(eqs)) {
      if ((mean.h2 - mean.h1)[eqs[i]] != 0) {
        rh[eqs[i]] <- as.numeric(mean.h1 < mean.h2)[eqs[i]] + 1
        uneq <- c(uneq, i)
      }
      if (length(uneq) > 0) { eqs <- eqs[-uneq] }
    }
    if (length(eqs) > 0) {
      warning(paste("Means are also equal for ",
                    case.ids[eqs],
                    "; choosing haplotype 1 by default", sep = ""))
      rh[eqs] <- 1
    }
  }

  risk.haps <- paste(case.ids, rh, sep = ".")
  return(risk.haps)
}

