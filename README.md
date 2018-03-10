# chavi
Consensus haplotype detection and visualization 

The R package chavi identifies and visualizes long shared IBS between case haplotypes. 
It is intended to be used in the setting where a particular genomic site is of interest 
(e.g., a mutation site) and the region around that site is expected to be highly conserved among carriers. 
We also assume that only one haplotype per case subject carries the feature. 

See https://aplantin.github.io/chavi/ for further documentation and sample usage. 

# Installation instructions: 

  % library(devtools)  <br/>
  install_github("aplantin/chavi")  <br/>
  library(chavi)  <br/>
