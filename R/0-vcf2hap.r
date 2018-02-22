#' vcf2hap
#'
#' This function loads a phased vcf file and a case ID file and produes
#' a haplotype file in required format.
#'
#' The vcf file should have subject identifiers (column names) that match the
#' case ID file, its second column should contain the genomic position (in BP
#' or cM, to match the position of the start and end sites), and the header
#' row should begin with '#CHROM'. This file may contain extra subjects; the
#' hap file produced will contain only the cases.
#'
#' The case ID file should be a text file with a single column containing all
#' of the case subject identifiers.
#'
#' The resulting haplotype file has CaseID1.1 and CaseID1.2 as the haplotype
#' identifiers for a subject with identifier CaseID1.
#'
#' @importFrom utils file_test read.table
#'
#' @param vcf.file Path to the input vcf file
#' @param case.id.file Path to the input case ID file
#' @return A matrix of case haplotypes
#' @export
#'
vcf2hap <- function(vcf.file, case.id.file){

  ## phased vcf, header not read in here
  vcf <- read.table(vcf.file, as.is = TRUE)

  ## get header by reading comment lines one at a time
  con <- file(vcf.file, "r")
  headline = FALSE
  while (headline == FALSE) {
    lin <- readLines(con, n = 1)
    bit <- substr(lin, start = 1, stop = 6)
    headline <- (bit == "#CHROM")
  }
  close(con)

  ## add header (via column names) to vcf
  colnames(vcf) <- strsplit(lin, "\t")[[1]]

  ## restrict to cases
  case.ids <- read.table(case.id.file, as.is = TRUE)$V1
  vcf <- vcf[, c(2, which(colnames(vcf) %in% case.ids))]

  ## split into haplotypes
  haps <- transform(vcf, lapply({
    l <- as.list(vcf[,-1], stringsAsFactors = FALSE)
    names(l) <- colnames(vcf)[-1]
    l
  }, function(x) {
    do.call(rbind, strsplit(x, '|', fixed  = TRUE)) }),
  stringsAsFactors = FALSE)
  haps <- haps[, -c(2:ncol(vcf))]

  ## return haplotype file
  return(haps)
}

