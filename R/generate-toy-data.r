#' Generates simulated data for testing CHap functions.
#'
#' @importFrom stats runif
#' @importFrom utils write.table
#'
#' @param seed Random seed for data generation; default 1.
#' @param n.markers Number of markers on simulated "chromosome"; default 50.
#' @param n.subj Number of subjects; default 10, all case subjects.
#' @param marker.spacing Distance between markers; default 0.1 (as if cM positions).
#' @param num.splits Number of times the consensus haplotype splits on each side
#' @keywords internal
#' @return Produces 3 files in the current directory: sim-case-ids.txt, sim-vcf.txt, and sim-ibs.txt. Returns R objects "haps" (the haplotype data frame for input into consensus) and "riskhaps" (the vector of true risk haplotypes)
#' @export
#'
generate_toy_data <- function(seed = 1000, n.markers = 50, n.subj = 10,
                              marker.spacing = 0.1, num.splits = 5) {
  set.seed(seed)

  ## Site of interest (this occurs between markers, but is allowed to cover a marker)
  site.start <- n.markers*marker.spacing/2 - 0.7*marker.spacing
  site.end <- n.markers*marker.spacing/2 - 0.3*marker.spacing

  ## Case ID file: All subjects are cases
  case.ids <- paste("SUBJID", 1:n.subj, sep = "")
  write.table(case.ids, file = "./sim-case-ids.txt", quote = FALSE,
              sep = "\n", row.names = FALSE, col.names = FALSE)


  ## Baseline haplotype file
  hapfile <- cbind(marker.spacing * 1:n.markers,
                   matrix(sample(c(0,1), size = n.markers * n.subj * 2, replace = T),
                          nrow = n.markers, ncol = n.subj*2))
  colnames(hapfile) <- c("POS", sapply(case.ids, FUN = function(x) paste(x, 1:2, sep = ".")))


  ## Choose which will be risk haplotype
  riskhaps <- sample(c(1:2), size = n.subj, replace = TRUE)


  ## Baseline shared risk hap
  true.consensus <- list()
  true.consensus[[1]] <- sample(c(0,1), size = n.markers, replace = TRUE)
  leftpos <- max(which(hapfile[,1] < site.start))
  rightpos <- min(which(hapfile[,1] > site.end))
  for (i in 1:num.splits) {
    if (i %% 2 == 0) {
      new.left <- floor(runif(1, 1, leftpos))
      new.right <- ceiling(runif(1, rightpos, n.markers))
      next.consensus <- true.consensus[[(i-1)]]
      next.consensus[1:new.left] <- sample(c(0,1), size = new.left, replace = TRUE)
      next.consensus[new.right:n.markers] <- sample(c(0,1), replace = TRUE,
                                                    size = (n.markers - new.right + 1))
      true.consensus[[(i+1)]] <- next.consensus
    } else {
      new.left <- floor(runif(1, 1, leftpos))
      new.right <- ceiling(runif(1, rightpos, n.markers))
      next.consensus <- true.consensus[[i]]
      next.consensus[1:new.left] <- sample(c(0,1), size = new.left, replace = TRUE)
      next.consensus[new.right:n.markers] <- sample(c(0,1), replace = TRUE,
                                                    size = (n.markers - new.right + 1))
      true.consensus[[(i+1)]] <- next.consensus
    }
    leftpos = new.left
    rightpos = new.right
  }


  ## Give each subject a risk haplotype (to right and left of site)
  leftpos <- max(which(hapfile[,1] < site.start))
  rightpos <- min(which(hapfile[,1] > site.end))
  hapsL <- sample(rep(1:num.splits, each = (num.splits %% n.subj + 1)))[1:n.subj]
  hapsR <- sample(rep(1:num.splits, each = (num.splits %% n.subj + 1)))[1:n.subj]
  for (i in 1:n.subj) {
    this.rh <- riskhaps[[i]]
    this.endL <- floor(runif(1, 1, leftpos))
    this.endR <- ceiling(runif(1, rightpos, n.markers))
    hapL <- hapsL[i] + 1
    hapR <- hapsR[i] + 1

    this.riskhap <- paste("SUBJID", i, ".", this.rh, sep = "")
    hapseq <- hapfile[, this.riskhap]
    hapseq[ this.endL:leftpos ] <- true.consensus[[hapL]][ this.endL:leftpos ]
    hapseq[ rightpos:this.endR ] <- true.consensus[[hapR]][ rightpos:this.endR ]
    hapfile[, this.riskhap] <- hapseq
  }


  ## Generate IBS file
  ibs <- data.frame(id1=character(),
                    hap1=numeric(),
                    id2=character(),
                    hap2=numeric(),
                    chrom=integer(),
                    start=numeric(),
                    end=numeric(),
                    length=numeric(),
                    stringsAsFactors=FALSE)

  for (i in 1:(2*n.subj-1)) {
    for (j in (i+1):(2*n.subj)) {
      id1 = ceiling(i/(2*n.subj)*10)
      hap1 = 2 - i%%2
      id2 = ceiling(j/(2*n.subj)*10)
      hap2 = 2 - j%%2
      chrom = 1
      hapseq1 = hapfile[, paste("SUBJID", id1, ".", hap1, sep = "")]
      hapseq2 = hapfile[, paste("SUBJID", id2, ".", hap2, sep = "")]
      same <- 1 - abs(hapseq1 - hapseq2)   ## 1 is same
      # find IBS segments
      startpos <- min(which(same == 1))
      endpos <- startpos
      while (endpos < n.markers) {
        if (any(same[startpos:n.markers] == 0)) {  # hap ends before "chromosome" ends
          endpos <- startpos + min(which(same[startpos:n.markers] == 0)) - 1
          ibs <- rbind(ibs, data.frame(id1 = case.ids[id1], hap1,
                                       id2 = case.ids[id2], hap2,
                                       chrom,
                                       start = hapfile[startpos,1],
                                       end = hapfile[endpos,1],
                                       length = NA))
          if (any(same[endpos:n.markers] == 1)) {
            startpos <- endpos + min(which(same[endpos:n.markers] == 1)) - 1  # start next segment
          } else {
            endpos <- n.markers  # no more segments
          }
        } else {
          endpos <- n.markers
          ibs <- rbind(ibs, data.frame(id1 = case.ids[id1], hap1,
                                       id2 = case.ids[id2], hap2,
                                       chrom,
                                       start = hapfile[startpos,1],
                                       end = hapfile[endpos,1],
                                       length = NA))
        }
      }
    }
  }
  ibs$length <- round(ibs$end - ibs$start, 1)
  keep <- which(ibs$length > 0.1)
  long.ibs <- ibs[keep, ]
  rownames(long.ibs) <- NULL
  write.table(long.ibs, file = "./sim-ibs.txt", sep = "\t", quote = F, row.names = F, col.names = F)


  ## Turn hap file into a VCF file (for demonstration)
  vcfhaps <- data.frame(sapply(1:n.subj, FUN = function(x) {
    paste(hapfile[, paste("SUBJID", x, ".1", sep = "")],
          hapfile[, paste("SUBJID", x, ".2", sep = "")],
          sep = "|") }))
  colnames(vcfhaps) <- case.ids

  toy.vcf <- data.frame(
    CHROM = rep(1,n.markers),
    POS = marker.spacing * c(1:n.markers),
    ID = paste("rs000", c(1:n.markers), sep = ""),
    REF = "G",
    ALT = "A",
    QUAL = ".",
    FILTER = "PASS",
    FORMAT = ".",
    vcfhaps,
    stringsAsFactors = FALSE
  )

  con <- file("./sim-haps.vcf", open="wt")
  writeLines("##INFO=<This.is.a.simulated.vcf.file>", con)
  writeLines(paste(c("#CHROM", colnames(toy.vcf)[-(1)]), collapse = "\t"), con)
  write.table(toy.vcf, con, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  close(con)

  haps <- hapfile
  riskhaps <- paste(case.ids, riskhaps, sep = ".")

  return(list(haps = haps, riskhaps = riskhaps))
}

