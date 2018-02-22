#' plot_consensus
#'
#' Produces a plot of consensus haplotypes. Haplotype segments belonging to a
#' single group are identified by color and (optionally) symbol. Site of interest
#' may be plotted; start and end locations appear as solid lines. Genomic
#' positions where consensus haplotype splits into subgroups may also be plotted;
#' these appear as dotted lines. Subjects may be ordered by group membership with
#' group names included alongside the subject identifiers.
#'
#' @importFrom graphics abline axis par plot points rect
#'
#' @param consensus Data frame containing consensus haplotype segments (columns are riskhap ID, segment start, segment end, and consensus haplotype segment ID. Output from get.consensus().
#' @param splits Vector of genomic positions where consensus haplotype splits into groups. Default NULL; if included, dashed lines will be drawn at each split position. Output frum get.consensus().
#' @param site.start Genomic position of start of site of interest. Default NULL; if included, a solid line will be drawn at site.start.
#' @param site.end Genomic position of end of site of interest. Default NULL; if included, a solid line will be drawn at site.end.
#' @param colors Vector of colors for haplotype segments. Default NULL (program will generate colors).
#' @param ncolors Number of unique colors to include on plot. Max 9, default 9. Colors are generated from brewer.pal to be as easily distinguishable as possible and repeated if more haplotype segments are present than colors.
#' @param symbols Logical indicating whether symbols should be overlaid line segments for each consensus haplotype segment. Default TRUE.
#' @param groups Vector identifying group membership for each subject. Default NULL.
#' @param group.names Names of subgroups of subjects. Default NULL.
#' @param label.size Controls relative size of subject and group labels. Default 1.
#' @param title.size Controls relative size of text in title and axis labels. Default 1.
#' @return Creates plot
#' @export
#'
plot_consensus <- function(consensus,
                           splits = NULL,
                           site.start = NULL, site.end = NULL,
                           colors = NULL, ncolors = 9,
                           symbols = TRUE,
                           groups = NULL, group.names = NULL,
                           label.size = 1, title.size = 1){
  if (!is.null(groups)) {
    consensus$id <- factor(consensus$id,
                           levels = consensus$id[order(as.numeric(as.factor(groups)))])
    if (is.null(group.names)) {
      group.names = levels(groups)
    }
  }
  consensus$id <- as.factor(consensus$id)  # need ID to be factor variable


  ## find plot limits
  limits <- range(c(consensus$start, consensus$end))
  rg <- limits[2] - limits[1]
  limits <- limits + c(-1,1) * rg/25
  limits[1] <- max(0, limits[1])

  ## set up color palette
  p <- length(unique(consensus$hapID))
  if (is.null(colors)) {
    mycols <- rep(brewer.pal(ncolors, "Set1"), p*(p%%ncolors + 1))[1:p]
  } else {
    ncolors = length(colors)
    mycols = rep(colors, p*(p%%ncolors + 1))[1:p]
  }

  ## Symbols
  mysymb <- rep(c(1:19), p*(p%%19 + 1))[1:p]

  ## create plot
  id.width <- 0.67*max(nchar(unique(as.character(consensus$id))))
  id.width = id.width*label.size
  if (!is.null(groups)) { id.width = id.width + 2 }
  par(mar = c(5.1, 2.1, 4.1, id.width+0.1))
  plot(rep(0, 2) ~ limits, type='n', xlab="Position", yaxt='n',
       ylab="", ylim=c(0.5, length(unique(consensus$id)) + 0.5),
       main="Consensus Risk Haplotypes",
       cex.main=title.size, cex.lab=title.size)
  for(i in 1:nrow(consensus)){
    this <- consensus[i,]
    rect(xleft=this$start, ybottom=as.numeric(this$id)-0.05,
         xright=this$end, ytop=as.numeric(this$id)+0.05,
         col=mycols[this$hapID], border=NA)
    if (symbols) {
      xpos = 0.5*(this$end + this$start)
      ypos = as.numeric(this$id)
      points(x = xpos, y = ypos, pch = mysymb[this$hapID],
             cex = label.size, col = "black", lwd = 2) }
  }

  ## if desired, add vertical dotted lines for splits
  if(!is.null(splits)){
    for(i in 1:length(splits)){
      abline( v=splits[i], lty=2 )
    }
  }

  ## if desired, add vertical solid line at start and end of site of interest
  if (!is.null(site.start)) {
    abline(v = site.start, col = "black", lwd = 1, lty = 1)
  }
  if (!is.null(site.end)) {
    abline(v = site.end, col = "black", lwd = 1, lty = 1)
  }

  ## label groups (if included)
  if (!is.null(groups)) {
    group.sizes <- table(groups)
    group.starts <- 0.75 + cumsum(c(0, group.sizes[-length(group.sizes)]))
    group.ends <- 0.25 + cumsum(group.sizes)
    axis(4, at=(group.starts + group.ends)/2, labels = group.names,
         pos=limits[2], lty=0, padj=0, cex.axis=label.size)
    axis(4, at=as.numeric(consensus$id), labels = consensus$id,
         line = 1.5, lty=0, las=2, cex.axis=label.size)
    for (i in 1:length(group.sizes)) {
      axis(4, tick = T, lwd.ticks = 0, line = 1.5,
           labels = F, at = c(group.starts[i], group.ends[i]))
    }
  } else {
    ## label axis with IDs in consensus file
    axis(4, at=as.numeric(consensus$id), labels = consensus$id,
         pos=limits[2]+rg/25, lty=0, las=2, cex.axis=label.size)
  }
}



