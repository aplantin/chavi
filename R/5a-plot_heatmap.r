#' plot_heatmap
#'
#' Generates heatmap to visualize pairwise IBS sharing between case subjects
#' at the site of interest
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices rainbow
#' @importFrom gridExtra grid.arrange
#' @importFrom lattice levelplot panel.levelplot panel.abline
#'
#' @param consensus Data frame containing consensus haplotype segments (columns are riskhap ID, segment start, segment end, and consensus haplotype segment ID. Output from get.consensus().
#' @param groups Vector indicating group membership. Default NULL; if included, subjects will be ordered by group (i.e., subjects in same group will be plotted next to each other.)
#' @param group.lines Logical indicating whether groups should be marked on plot. Default FALSE; if TRUE, gray lines will demarcate groups.
#' @param graphmin Minimum IBS value for plot. Anything lower (but nonzero) will be increased to this value for the plot. Default NULL.
#' @param graphmax Maximum IBS value for plot. Anything higher willl be truncated to this value for the plot. Default NULL.
#' @param main Title of plot
#' @param cx Multiplier for font size for subject labels, graph title. Default 0.75.
#' @return Creates IBS heat map.
#' @export
#'
ibs_heatmap <- function(consensus,
                        groups = NULL, group.lines = FALSE,
                        graphmin = NULL, graphmax = NULL,
                        main = NULL, cx = 0.75){

  ibs.mat = ibs_matrix(consensus, groups)
  trunc.out = trunc_ibs_mat(ibs.mat, graphmin, graphmax)

  if (is.null(graphmin)) {
    graphmin <- min(trunc.out[trunc.out != 0])*0.95
  }
  if (is.null(graphmax)) {
    graphmax <- max(trunc.out)*1.05
  }

  mycol = rainbow(200, start=0.2, s=log(c(1,10:208))/log(208))

  if (!is.null(groups)) {
    group.len <- table(groups)
    group.sep = cumsum(group.len) + 0.5
  }

  plot1 = levelplot(trunc.out, panel = function(...) {
    panel.levelplot(...)
    if (!is.null(groups)) {
      for (i in 1:length(unique(groups))) {
        panel.abline(h = group.sep[i])
        panel.abline(v = group.sep[i])
      }
    }
  },
  main = list(label = main, cex = cx),
  col.regions = mycol,
  at = seq(graphmin, graphmax, length.out = 200),
  ylab = "", xlab = "",
  scales = list(x = list(rot = 90, cex = cx), y = list(cex = cx)),
  colorkey = list(space = "left", labels = list(cex = cx))
  )
  grid.arrange(plot1)
}


