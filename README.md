# chavi
Consensus haplotype detection and visualization 

chavi identifies and visualizes long shared identity-by-state (IBS) segments 
      between case haplotypes. It is intended to be used in the setting where a 
      particular genomic site is of interest (e.g., a mutation site) and the region 
      around that site is expected to be highly conserved among carriers. The 
      similarity and length of this shared segment may be used to better understand 
      the relationships among carriers and the age of the genomic feature.

# Installation instructions: 

  % library(devtools)  <br/>
  install_github("aplantin/chavi")  <br/>
  library(chavi)  <br/>
