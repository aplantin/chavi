#' ibs_matrix
#'
#' Generates matrix of pairwise IBS sharing between case subjects
#' at the site of interest
#'
#' @param consensus Data frame containing consensus haplotype segments (columns are riskhap ID, segment start, segment end, and consensus haplotype segment ID. Output from get.consensus().
#' @param groups Vector indicating group membership.
#' @return Matrix of pairwise IBS sharing between case individuals across site of interest.
#' @export
#'
ibs_matrix <- function(consensus, groups = NULL) {
  ## Matrix of shared IBD lengths
  riskhaps <- unique(consensus$id)
  if (!is.null(groups)) {
    if (length(groups) != length(riskhaps)) {
      stop("Not every case subject is assigned a group!")
    }
    riskhaps <- riskhaps[order(as.numeric(as.factor(groups)))]
  }
  nhaps <- length(riskhaps)
  ibs.mat <- matrix(nrow = nhaps, ncol = nhaps)
  rownames(ibs.mat) = colnames(ibs.mat) = riskhaps

  for (i in 1:nhaps) {
    for (j in i:nhaps) {
      if (i == j) {
        ibs.mat[i, j] <- 0
      } else {
        ## longest pairwise hap for i and j
        res <- pairwise_summary(consensus,
                                as.character(riskhaps[i]),
                                as.character(riskhaps[j]))

        ## update in matrix
        ibs.mat[i, j] = ibs.mat[j, i] = res$length
      }
    }
  }
  return(ibs.mat)
}
