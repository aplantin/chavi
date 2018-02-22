## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  fig.width = 7, fig.height = 3.5, 
  comment = "#>"
)
library(devtools)
install_github("aplantin/chavi")
library(chavi)

## ----load-data-----------------------------------------------------------
vcf.file <- system.file("extdata", "sim-haps.vcf", package = "chavi")
case.id.file <- system.file("extdata", "sim-case-ids.txt", package = "chavi")
haps <- vcf2hap(vcf.file, case.id.file)

## ----view-haps-----------------------------------------------------------
haps[1:5,1:5]

## ----find-risk-----------------------------------------------------------
ibs.file <- system.file("extdata", "sim-ibs.txt", package = "chavi")
riskhaps <- id_riskhap(ibs.file, case.id.file, site.start = 2.43, site.end = 2.47)
riskhaps

## ----consensus-----------------------------------------------------------
cons.out <- get_consensus(haps, riskhaps, site.start = 2.43, site.end = 2.47, verbose = FALSE)

## ----consensus-1---------------------------------------------------------
head(cons.out$consensus)

## ----consensus-summ------------------------------------------------------
consensus_summary(cons.out$consensus)

## ----consensus-2---------------------------------------------------------
head(cons.out$splits)

## ----plot-consensus-1----------------------------------------------------
plot_consensus(cons.out$consensus, symbols = F, 
               label.size = 0.75)

## ----plot-consensus-2----------------------------------------------------
plot_consensus(cons.out$consensus, symbols = T, 
               label.size = 0.75)

## ----plot-consensus-3----------------------------------------------------
plot_consensus(cons.out$consensus, site.start = 2.43, site.end = 2.47, symbols = T, 
               label.size = 0.75)

## ----plot-consensus-4----------------------------------------------------
plot_consensus(cons.out$consensus, cons.out$splits, 
               site.start = 2.43, site.end = 2.47, symbols = T, 
               label.size = 0.75)

## ----plot-consensus-5----------------------------------------------------
plot_consensus(cons.out$consensus, cons.out$splits, 
               site.start = 2.43, site.end = 2.47, symbols = T, 
               colors = c("red", "blue", "green"), 
               label.size = 0.75)

## ----plot-consensus-6----------------------------------------------------
plot_consensus(cons.out$consensus, cons.out$splits, 
               site.start = 2.43, site.end = 2.47, symbols = T, 
               groups = c(2, 1, 1, 2, 1, 2, 2, 1, 2, 2), 
               group.names = c("Pop 1", "Pop 2"), 
               label.size = 0.75)

## ----plot-heatmap, fig.width = 5, fig.height = 3.5-----------------------
ibs_heatmap(cons.out$consensus)

## ----plot-heatmap-2, fig.width = 5, fig.height = 3.5---------------------
ibs_heatmap(cons.out$consensus, groups = c(2, 1, 1, 2, 1, 2, 2, 1, 2, 2), 
            group.lines = T)

