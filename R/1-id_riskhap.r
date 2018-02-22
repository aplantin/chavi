#' id_riskhap
#'
#' Identify which haplotype is the risk haplotype for each case individual.
#' Note that we assume only one haplotype carries the feature at the site
#' of interest.
#'
#' This function loads an IBS file and a case ID file, and also requires the
#' start and end locations of the site of interest to be specified.
#'
#' The IBS file should have eight columns: ID1, Hap1, ID2, Hap2, Chrom,
#' Start, End, Length. It may be generated using Brian Browning's ibs.jar
#' utility, located at [].
#'
#' The case ID file should be a text file with a single column containing all
#' of the case subject identifiers.
#'
#' site.start and site.end are the genomic positions that begin and end the
#' site of interest, in the same units as those used in the IBS file and the
#' vcf/haplotype file.
#'
#' Of the two haplotypes for a given case, the one with the longest median IBS
#' length at the site of interest (shared with other case subjects) is taken
#' to be the "risk haplotype" and carry the feature of interest. If the median
#' length is equal for the two haplotypes, the mean is used as a tie-breaker. If
#' the means are also equal, then a warning is generated and by default
#' haplotype 1 is identified as the risk haplotype.
#'
#' @importFrom stats aggregate median
#' @importFrom utils file_test read.table
#' @param ibs.file Path to the input IBS file
#' @param case.id.file Path to the input case ID file
#' @param site.start Genomic position at start of site of interest
#' @param site.end Genomic position at end of site of interest
#' @return A vector of risk haplotype IDs
#' @export
#'
id_riskhap <- function(ibs.file, case.id.file, site.start, site.end){

  ## read and check IBS
  if (file_test("-f", ibs.file)) {
    ibs <- read.table(ibs.file, as.is=TRUE, header=FALSE)
  } else {
    stop("IBS file not found.")
  }
  if (ncol(ibs) != 8) {
    stop("Please input IBS file with eight columns: ID.1, HapID.1, ID.2,
         HapID.2, Chromosome, Start, End, Length")
  }
  colnames(ibs) <- c("id1", "hap1", "id2", "hap2", "chrom", "start", "end",
                     "length")

  ## read and check case IDs
  if (file_test("-f", case.id.file)) {
    case.ids <- read.table(case.id.file, as.is = TRUE)$V1
  } else {
    stop("Case ID file not found.")
  }
  if (length(case.ids) == 0) {
    stop("No case IDs found. Please make sure the case ID file is a text file
         containing a single column of case identifiers.")
  }
  if (all(!case.ids %in% c(ibs$id1, ibs$id2))) {
    stop("No IBS was found for the case IDs specified. Please check that the
         case identifiers match the subject IDs in the IBS file.")
  } else if (any(!case.ids %in% c(ibs$id1, ibs$id2))) {
    warning("Some case IDs not found in the IBS file. If IBS is expected for all
            cases, please check that the case identifiers match the subject IDs
            in the IBS file.")
  }

  ## Check that site.start and site.end are reasonable (i.e., same units as IBS)
  if (max(ibs$end) > 300 & site.end < 300) {  # 300 > cM length of chromosome 1
    warning(paste("The IBS start and end locations range from ", min(ibs$start),
                  " to ", max(ibs$end), ", and the site of interest is (",
                  site.start, ", ", site.end, "). Please double check that your
                  IBS locations and site position are in the same units.",
                  sep = ""))
  }

  out <- id_riskhap_inner(ibs, case.ids, site.start, site.end)

  return(out)
}


