# chavi

The R package chavi performs consensus haplotype detection and visualization. See vignette for sample usage on simulated data. 

# Installation instructions: 

  library(devtools)  <br/>
  install_github("aplantin/chavi")  <br/>
  library(chavi)  <br/>


```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  fig.width = 7, fig.height = 3.5, 
  comment = "#>"
)
library(devtools)
install_github("aplantin/chavi")
library(chavi)
```

# Description 

chavi is a package that identifies and visualizes long shared IBS between case haplotypes. It is intended to be used in the setting where a particular genomic site is of interest (e.g., a mutation site) and the region around that site is expected to be highly conserved among carriers. We also assume that only one haplotype per case subject carries the feature. 

Vocabulary used throughout this document: 

- "Case subjects" are those that carry the site of interest on one chromosome copy. In the case where the feature of interest is a mutation known to be causal for a disease, case subjects are those with the disease. 
- "Consensus risk haplotype" refers to the stretch of genome shared identically between (nearly) all cases around the site of interest. 
- "Risk haplotype" refers to the copy of the chromosome for each subject that carries the feature of interest. 

## Step 1: Load and prepare data 

We assume that the available data include a phased VCF file and a text file with one case subject identifier per line. The VCF file may include other subjects in addition to cases; these will be filtered out for the remainder of the analysis. 

The VCF file is converted to a "haplotype" file. Rows are genetic markers, and there is a column labeled POS identifying the position of the genetic markers (in bp or cM). All of the markers must be on a single chromosome. The remaining columns contain the data for each case haplotype. The haplotype identifiers are SUBJ.1 and SUBJ.2 for the first and second copy of the chromosome for the case subject with identifier SUBJ. Other columns may be included and will be ignored.  


```{r load-data}
vcf.file <- system.file("extdata", "sim-haps.vcf", package = "chavi")
case.id.file <- system.file("extdata", "sim-case-ids.txt", package = "chavi")
haps <- vcf2hap(vcf.file, case.id.file)
```

The haplotype file looks like this: 
```{r view-haps}
haps[1:5,1:5]
```


After the haplotype file is loaded, we need to identify which copy of the chromosome carries the feature of interest for each case subject. This procedure assumes there is long IBS shared among case subjects at that site. The data were generated such that the site of interest is from 2.43 to 2.47 cM. 

```{r find-risk} 
ibs.file <- system.file("extdata", "sim-ibs.txt", package = "chavi")
riskhaps <- id_riskhap(ibs.file, case.id.file, site.start = 2.43, site.end = 2.47)
riskhaps
```

This function tells us that Haplotype 2 of Subjects 1, 2, 3, and 6 carry the feature of interest, whereas for the remaining subjects, Haplotype 1 carries the feature. 


## Step 2: Identify consensus risk haplotype


```{r consensus}
cons.out <- get_consensus(haps, riskhaps, site.start = 2.43, site.end = 2.47, verbose = FALSE)
``` 


The consensus haplotype file contains [...]. 
```{r consensus-1}
head(cons.out$consensus)
```

We can summarize this data frame as the longest haplotype carried by each case subject: 
```{r consensus-summ}
consensus_summary(cons.out$consensus)
``` 

The split locations are also stored in the output of get_consensus(). 
```{r consensus-2}
head(cons.out$splits)
```


## Step 3: Visualize consensus risk haplotype 

The plot below demonstrates most of the plot options available for visualizing the consensus risk haplotypes. 

```{r plot-consensus-6}
plot_consensus(cons.out$consensus, cons.out$splits, 
               site.start = 2.43, site.end = 2.47, symbols = T, 
               groups = c(2, 1, 1, 2, 1, 2, 2, 1, 2, 2), 
               group.names = c("Pop 1", "Pop 2"), 
               label.size = 0.75)
```



An alternative (or complementary) visualization is a heat map of pairwise IBS segment lengths. The values on the heat map scale are segment lengths, so red colors indicate long pairwise segments and blue colors indicate shorter segments. Again, subjects may be organized into groups. 

```{r plot-heatmap-2, fig.width = 5, fig.height = 3.5} 
ibs_heatmap(cons.out$consensus, groups = c(2, 1, 1, 2, 1, 2, 2, 1, 2, 2), 
            group.lines = T)
```

